#' Open a DuckDB connection with views over pool/pointer parquets
#'
#' Creates two SQL views (`pointers`, `pool`) on a directory of per-sample
#' parquets produced by [Ms2PseudoRaw()]. Caller is responsible for
#' `DBI::dbDisconnect(con, shutdown = TRUE)` — or pass the connection to
#' [clean_pseudo_spectra()] which closes it on exit.
#'
#' @param file_directory Directory containing `*_pool.parquet` + `*_pointers.parquet`.
#' @return A live DuckDB connection.
#' @keywords internal
#' @export
setup_connection <- function(file_directory) {
  # Resolve to absolute path — DuckDB may not share R's working directory
  abs_dir <- normalizePath(file_directory, mustWork = TRUE)
  # Fail with something readable. DuckDB's own message for an empty glob is an
  # opaque "IO Error: No files found that match the pattern ...", which hides
  # the real cause: no sample produced pointers (e.g. every target fell outside
  # the RT window, or the mzML was truncated).
  if (!length(list.files(abs_dir, pattern = "_pointers\\.parquet$")))
    cli::cli_abort(c(
      "No {.file *_pointers.parquet} in {.path {abs_dir}}.",
      "i" = "No sample produced pointer rows - check that the targets' RT
             values fall inside the acquisition and that the mzML files are
             complete."
    ))
  # Use a local temp dir for DuckDB spill files to avoid /tmp exhaustion on shared nodes
  tmp_dir <- file.path(abs_dir, ".duckdb_tmp")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  con <- DBI::dbConnect(duckdb::duckdb(), config = list(temp_directory = tmp_dir))
  DBI::dbExecute(con, sprintf(
    "CREATE OR REPLACE VIEW pointers AS SELECT * FROM read_parquet('%s/*_pointers.parquet')", 
    abs_dir
  ))
  DBI::dbExecute(con, sprintf(
    "CREATE OR REPLACE VIEW pool AS SELECT * FROM read_parquet('%s/*_pool.parquet')", 
    abs_dir
  ))
  con
}

#' Fetch raw MS2 + MS1 rows for a vector of features from DuckDB
#'
#' Joins `pointers ⨝ pool` on `(sample, scan + 1)` with the precursor cut
#' `m.mz <= p.mz_target + 2`. Chunks the feature list to bound `IN (…)`.
#'
#' @param feature_ids Character vector of feature ids.
#' @param con A connection from [setup_connection()].
#' @param chunk_size Integer. Features per `IN`-clause.
#' @return `data.table` with columns
#'   `sample, feature, mz, fragment_mz, int, ms1_scan, ms2_scan, ms1_int`.
#' @keywords internal
#' @export
batch_get_feature_data <- function(feature_ids, con, chunk_size = 500) {
  # DuckDB is fast at IN clauses but we chunk to avoid overly long SQL
  chunks <- split(feature_ids, ceiling(seq_along(feature_ids) / chunk_size))
  
  results <- data.table::rbindlist(lapply(chunks, function(ids) {
    id_list <- paste0("'", ids, "'", collapse = ", ")
    query <- sprintf("
      WITH feat_pointers AS (
        SELECT * FROM pointers WHERE feature IN (%s)
      )
      SELECT
        p.sample,
        p.feature,
        p.mz_target AS mz,
        m.mz AS fragment_mz,
        m.int,
        p.scan AS ms1_scan,
        m.scan AS ms2_scan,
        p.ms1_int
      FROM feat_pointers p
      LEFT JOIN pool m ON (
        p.scan + 1 = m.scan AND
        p.sample = m.sample
      )
      WHERE (m.mz <= p.mz_target + 2 OR m.mz IS NULL)
      ORDER BY p.feature, p.sample, p.scan
    ", id_list)
    data.table::as.data.table(DBI::dbGetQuery(con, query))
  }))
  
  return(results)
}



#' Consensus m/z for one fragment group
#'
#' Intensity-weighted mean over the FWHM core — the centroids at or above
#' half the group's maximum intensity. Low-intensity centroids carry by far
#' the largest m/z error (per-point precision scales as `1/sqrt(I)`), so
#' they are dropped rather than averaged in. Falls back to the plain mean of
#' the core when the weights vanish, and to `NA_real_` when the group holds
#' no finite values.
#'
#' @keywords internal
#' @noRd
.consensus_mz <- function(mz, int, Isat = Inf, frac = 0.5) {
  ok <- is.finite(mz) & is.finite(int)
  if (!any(ok)) return(NA_real_)
  mz <- mz[ok]; int <- int[ok]
  keep <- int >= frac * max(int)     # FWHM core: centroids above half-max
  w    <- pmin(int[keep], Isat)      # weight; clip at Isat to not over-weight a saturated apex
  if (sum(w) <= 0) return(mean(mz[keep]))
  sum(w * mz[keep]) / sum(w)
}

#' Build a feature's per-group statistics table from raw MS1+MS2 data.table
#'
#' Stages, in order:
#'
#' \enumerate{
#'   \item per-sample apex window (`buffer`);
#'   \item pass 1 (within sample): `grouper` clusters raw `fragment_mz`
#'     scoped to each sample separately, so this pass only ever sees
#'     within-sample scan-to-scan jitter, never between-sample
#'     calibration drift;
#'   \item collapse each `(sample, pass-1 group)` to one evidence-weighted
#'     point - intensity-weighted consensus m/z (via [.consensus_mz()])
#'     and *summed* intensity, so a sample with more corroborating scans
#'     contributes proportionally more evidence, not more competing
#'     points;
#'   \item pass 2 (across samples): `grouper` clusters those per-sample
#'     points, resolving between-sample calibration drift - a separate,
#'     larger-scale effect pass 1 never attempts to address. The result
#'     is joined back onto every original raw row; `int`/`ms1_int` on
#'     `dt` are never touched by either pass, so MS1/MS2 correspondence
#'     stays at full `(sample, ms1_scan)` resolution throughout;
#'   \item consensus m/z — one value per final group, written to
#'     `group_mz` by [.consensus_mz()]. From here on m/z is a constant
#'     label;
#'   \item statistics closure, which reduces each group to one row and
#'     never touches m/z again — it just carries `group_mz` through as
#'     `fragment_mz`.
#' }
#'
#' No thresholding — see [clean_pseudo_spectra()] for the composition with a
#' [filter-contract] closure.
#'
#' The `(group, sample, ms1_scan)` aggregation computed inside the
#' statistics closure is attached to the return value as attribute
#' `"grouped_signal"` so downstream code (diagnostic plots) can reuse it
#' without a second groupby.
#'
#' @param dt Output of [batch_get_feature_data()] restricted to one feature.
#' @param mz_target Numeric. Precursor m/z (drives the ppm tolerance).
#' @param grouper A grouper closure - see [grouper-contract]. Used for
#'   both passes.
#' @param group_statistics_fn A statistics closure - see [stats-contract].
#' @param buffer Integer. Scan window around each sample's MS1 apex.
#' @param consensus_mz_fn A consensus-mz closure - see
#'   [consensus-mz-contract]. `NULL` (default) uses
#'   [consensus_mz_weighted()] (`frac = 0.5, Isat = 1e4`) - the package's
#'   long-standing default behavior.
#' @return A per-group summary `data.table` (see [stats-contract]) with
#'   one row per group, plus a `grouped_signal` attribute
#'   (aggregated `(group, sample, ms1_scan)`) used by the diagnostic report.
#' @keywords internal
#' @export
process_pseudo_spec <- function(
  dt,
  mz_target,
  grouper             = ppm_upgma(ppm = 30),
  group_statistics_fn = correlation_statistics(),
  buffer              = 50L,
  consensus_mz_fn     = NULL
) {
  if (is.null(consensus_mz_fn)) consensus_mz_fn <- consensus_mz_weighted()
  if (nrow(dt) == 0) return(.empty_correlations())
  dt <- dt[!is.na(int) & !is.na(fragment_mz)]
  if (nrow(dt) == 0) return(.empty_correlations())

  dt[, apex := ms1_scan[which.max(ms1_int)], by = sample]
  dt <- dt[abs(ms1_scan - apex) < buffer]
  if (nrow(dt) == 0) return(.empty_correlations())

  # Pass 1 (within sample): cluster raw fragment_mz per sample - within-
  # sample scan-to-scan jitter is a different, smaller-scale noise source
  # than between-sample calibration drift, so this pass only sees the
  # former.
  dt[, group1_local := grouper(mz        = fragment_mz,
                               intensity = int,
                               target_mz = mz_target,
                               scan      = ms2_scan,
                               ms1_int   = ms1_int),
     by = sample]

  # Local ids are only unique WITHIN one sample's own grouper call (two
  # samples' "1" are unrelated fragments). .GRP over (sample, group1_local)
  # assigns one collision-free, table-wide-unique integer per distinct
  # combination, independent of row order - it is not a stable id across
  # separate calls, only a disposable join key for this one call.
  dt[, group1 := .GRP, by = .(sample, group1_local)]

  # Disposable per-sample summary, pass-2 input only: one point per
  # (sample, group1) - intensity-weighted consensus mz (the same
  # estimator used for the final report) and SUMMED intensity, so a
  # sample with more corroborating scans contributes proportionally more
  # evidence to pass 2, not more competing points. dt itself (and its
  # ms1_scan-level resolution, which group_statistics_fn needs) is
  # untouched by this - only a new `group` label gets joined back below.
  pass2_input <- dt[, .(
    fragment_mz = consensus_mz_fn(fragment_mz, int),
    intensity   = sum(int)
  ), by = .(sample, group1)]

  # Pass 2 (across samples): cluster the per-sample summary points, then
  # propagate the result back onto every original raw row via (sample,
  # group1).
  pass2_input[, group := grouper(mz        = fragment_mz,
                                 intensity = intensity,
                                 target_mz = mz_target)]
  dt[pass2_input, on = .(sample, group1), group := i.group]
  dt[, c("group1_local", "group1") := NULL]

  # consensus m/z — computed once per group, before any statistics. m/z is a
  # constant label from this point on.
  dt[, group_mz := consensus_mz_fn(fragment_mz, int), by = group]

  group_statistics_fn(dt)
}

#' Legacy per-group thresholder (superseded by [correlation_filtering()])
#'
#' Kept as an internal convenience for old scripts. New code should use
#' `correlation_filtering(cor_thresh, fdr_thresh)` — a [filter-contract]
#' closure that composes with [process_pseudo_spec()] via
#' [clean_pseudo_spectra()]'s `filter_fn` argument.
#'
#' @param correlations Per-group summary from a [stats-contract] closure
#'   (typically [correlation_statistics()]).
#' @param cor_thresh Numeric. Minimum correlation to keep.
#' @param fdr_thresh Numeric. Maximum FDR to keep.
#' @return A `data.table` (possibly empty).
#' @keywords internal
select_pseudo_spec <- function(correlations, cor_thresh = 0.5, fdr_thresh = 0.05) {
  if (nrow(correlations) == 0) return(correlations)
  correlations[correlation >= cor_thresh & fdr <= fdr_thresh]
}

#' Internal: diagnostic plots for one pseudo-spectrum (ggplot2-gated).
#' Returns `list()` when ggplot2 is unavailable or there is nothing to plot.
#'
#' `stats` is the full per-group stats table (before filtering); `pseudo_spec`
#' is the post-filter subset. `grouped_signal` is the per-`(group, sample,
#' ms1_scan)` aggregation typically pulled from `attr(stats, "grouped_signal")`.
#'
#' @keywords internal
#' @noRd
.pseudo_spec_plots <- function(grouped_signal, stats, pseudo_spec, feature_id) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(list())
  if (nrow(pseudo_spec) == 0) return(list())

  # Color contrast reference: high-int in the pre-filter distribution, so the
  # scale doesn't collapse when survivors happen to be uniformly bright.
  ref_int <- if (nrow(stats)) stats$int else pseudo_spec$int
  q95 <- stats::quantile(ref_int, 0.95, na.rm = TRUE)

  p_corr <- ggplot2::ggplot(pseudo_spec,
                            ggplot2::aes(x = correlation, y = fdr,
                                         color = int > q95)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::labs(title = paste0(feature_id, ": correlation vs FDR (surviving groups)"),
                  color = "high int") +
    ggplot2::theme_minimal(base_size = 9)

  p_spec <- ggplot2::ggplot(pseudo_spec) +
    ggplot2::geom_segment(ggplot2::aes(x = fragment_mz, xend = fragment_mz,
                                       y = 0, yend = int, color = 1 - fdr))
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_spec <- p_spec + ggrepel::geom_text_repel(
      data = pseudo_spec,
      ggplot2::aes(x = fragment_mz, y = int, label = round(fragment_mz, 3)),
      direction = "x", segment.size = 0.3, size = 2.5
    )
  }
  p_spec <- p_spec +
    ggplot2::labs(title = paste0(feature_id, ": pseudo-spectrum"),
                  x = "fragment m/z", y = "intensity") +
    ggplot2::theme_minimal(base_size = 9)

  filtered <- if (is.null(grouped_signal) || !nrow(grouped_signal)) NULL
              else grouped_signal[group %in% pseudo_spec$group]

  p_ms1 <- if (!is.null(filtered) && nrow(filtered)) {
    ms1_summary <- filtered[, .(ms1_int = mean(ms1_int)), by = .(sample, ms1_scan)]
    ggplot2::ggplot(ms1_summary, ggplot2::aes(x = ms1_scan, y = ms1_int)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ sample, ncol = 2, scales = "free") +
      ggplot2::labs(title = paste0(feature_id, ": MS1 chromatograms")) +
      ggplot2::theme_minimal(base_size = 8)
  } else NULL

  p_eic <- if (!is.null(filtered) && nrow(filtered)) {
    eic <- data.table::copy(filtered)[, int := scale(int)[, 1], by = .(sample, group)]
    ggplot2::ggplot(eic, ggplot2::aes(x = ms1_scan, y = int,
                                      color = factor(group))) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ sample, ncol = 2, scales = "free") +
      ggplot2::labs(title = paste0(feature_id, ": fragment EICs (scaled)")) +
      ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(legend.position = "none")
  } else NULL

  list(p_corr = p_corr, p_spec = p_spec, p_ms1 = p_ms1, p_eic = p_eic)
}

#' Build a single-row DeconvERROR sentinel matching the real-row schema
#'
#' Returned by [clean_pseudo_spectra()]'s per-target loop whenever a target
#' could not be turned into a real spectrum, so every requested target
#' contributes exactly one row (real fragment rows, or exactly one of
#' these) to the return value — never silently zero rows. `deconv_error`
#' is `NA` on real rows and the reason code here on sentinel rows;
#' [export_mgf()] reads it to decide whether to write a real MGF block or
#' a `DeconvERROR:<code>` placeholder for that feature.
#'
#' @keywords internal
#' @noRd
.deconv_error_row <- function(target, error_code, mz_target = NA_real_) {
  data.table::data.table(
    mz_g2           = NA_real_,
    grouped_mz      = NA_integer_,
    mz              = mz_target,
    total_intensity = NA_real_,
    n_samples       = NA_integer_,
    count           = NA_integer_,
    r               = NA_real_,
    t_stat          = NA_real_,
    p_value         = NA_real_,
    fdr             = NA_real_,
    feature_id      = target,
    deconv_error    = error_code
  )
}

#' Clean pseudo-spectra for many features (main batch entry)
#'
#' Orchestrates: connect → fetch → per-feature [process_pseudo_spec()]
#' (grouper + `group_statistics_fn`) → `filter_fn` → optional plot / MGF.
#' Returns a `data.table` in the schema downstream `03_make_mgf.R` expects.
#'
#' Guarantees exactly one outcome per element of `targets`: either real
#' fragment rows, or a single [.deconv_error_row()] sentinel (`deconv_error`
#' set to a reason code — `"MissingTargetInfo"`, `"MissingMs1Signal"`, or
#' `"NoAssociatedFragments"`). No target is ever silently dropped from the
#' return value; nothing here aborts on bad input, including a target
#' directory with no pointer parquets at all. This is what lets
#' [export_mgf()] guarantee one MGF block per `fd` row regardless of data
#' quality — the property [process_xcms()]/[process_targetlist()] batching
#' depends on, since each batch must independently account for its own
#' slice of targets.
#'
#' @param targets Character vector of feature ids.
#' @param file_directory Path to pool/pointer parquets.
#' @param fd Optional `featureDefinitions()`-shaped frame; supplies
#'   `mzmed` and `rtmed` per target. A target missing from `fd`, or whose
#'   `mzmed` (or, when `require_rt = TRUE`, `rtmed`) is `NA`, is reported
#'   as a `"MissingTargetInfo"` sentinel row rather than processed.
#' @param con Optional pre-opened connection (else opened + closed internally).
#' @param buffer,fdr_thresh,cor_thresh,min_count Tunables passed downstream.
#'   `min_count` is consumed by the default `group_statistics_fn`;
#'   `cor_thresh` / `fdr_thresh` are consumed by the default `filter_fn`.
#'   A custom `group_statistics_fn` / `filter_fn` may ignore them.
#' @param grouper A grouper closure - see [grouper-contract]. Used for
#'   both of [process_pseudo_spec()]'s passes (within-sample, then
#'   across-sample).
#' @param group_statistics_fn A statistics closure conforming to
#'   [stats-contract] — receives the raw per-feature dt (with `group`
#'   column) and returns the per-group summary (no thresholding).
#'   `NULL` (default) uses [correlation_statistics()] with `min_count`.
#'   Swap in [apex_correlation_statistics()] to add apex-shift columns.
#' @param consensus_mz_fn A consensus-mz closure - see
#'   [consensus-mz-contract]. `NULL` (default) uses
#'   [consensus_mz_weighted()] - the package's long-standing default.
#' @param filter_fn A filter closure conforming to [filter-contract] —
#'   receives the per-group stats table and returns a row-subset.
#'   `NULL` (default) uses [correlation_filtering()] with
#'   `cor_thresh`/`fdr_thresh`. Swap in [apex_correlation_filtering()]
#'   for apex-aware thresholding (pair with
#'   [apex_correlation_statistics()]).
#' @param require_rt Logical. When `TRUE`, a target whose `fd$rtmed` is
#'   `NA` is also reported as `"MissingTargetInfo"` (not just `mzmed`).
#'   Set by routes where RT drives scan extraction upstream (the
#'   target-list route); `FALSE` (default) for routes where RT is purely
#'   informational (the xcms route, where per-sample scan windows come
#'   from `chromPeaks`, not this feature-level summary).
#' @param batch_idx,n_batches Optional slice: process only batch
#'   `batch_idx` of `n_batches` equal contiguous chunks of `targets`.
#'   `batch_idx = NULL` (default) processes every target in-process.
#'   Set both for the DAG-batched `stage = "filter"` path.
#' @param report,report_path PDF diagnostic report on/off + path.
#' @param mgf_dir Optional directory for per-feature MGFs.
#' @param fetch_chunk_size Features per DuckDB query.
#' @param progress Print a per-feature progress counter.
#' @return A `data.table` matching the legacy 02_filter_features schema:
#'   `mz_g2, grouped_mz, mz, total_intensity, n_samples, count, r, t_stat,
#'   p_value, fdr, feature_id`, plus `deconv_error` (`NA` on real rows).
#'   `mz_g2` and `mz` both carry the per-group consensus m/z.
#'   `total_intensity` is a *mean* despite the name — the mean across
#'   `(sample, ms1_scan)` cells of each cell's summed fragment intensity,
#'   i.e. a representative peak height that does not scale with sample or
#'   scan count. Column names are fixed by the downstream
#'   `modules/annotation` contract; `deconv_error` is additive — consumers
#'   reading these rows directly (bypassing [export_mgf()]) should filter
#'   `is.na(deconv_error)` to get only real fragment rows.
#' @export
clean_pseudo_spectra <- function(
  targets,
  file_directory,
  fd                  = NULL,
  con                 = NULL,
  grouper             = ppm_upgma(ppm = 30),
  buffer              = 50L,
  fdr_thresh          = 0.01,
  cor_thresh          = 0.3,
  min_count           = 10L,
  group_statistics_fn = NULL,
  filter_fn           = NULL,
  consensus_mz_fn     = NULL,
  require_rt          = FALSE,
  batch_idx           = NULL,
  n_batches           = 1L,
  report              = TRUE,
  report_path         = "pseudo_specs_report.pdf",
  mgf_dir             = NULL,
  fetch_chunk_size    = 200L,
  progress            = TRUE
) {
  if (!is.null(batch_idx)) {
    per_batch <- ceiling(length(targets) / n_batches)
    my_start  <- 1L + (as.integer(batch_idx) - 1L) * per_batch
    my_end    <- min(as.integer(batch_idx) * per_batch, length(targets))
    if (my_start > length(targets)) return(data.table::data.table())
    targets <- targets[my_start:my_end]
  }

  if (is.null(group_statistics_fn))
    group_statistics_fn <- correlation_statistics(min_count = min_count)
  if (is.null(filter_fn))
    filter_fn <- correlation_filtering(cor_thresh = cor_thresh,
                                       fdr_thresh = fdr_thresh)

  fd_df <- NULL
  if (!is.null(fd)) {
    fd_df <- as.data.frame(fd)
    if ("feature_id" %in% colnames(fd_df)) rownames(fd_df) <- fd_df$feature_id
      if (!"mzmed" %in% colnames(fd_df) && "mz" %in% colnames(fd_df))
    fd_df$mzmed <- fd_df$mz
    if (!"rtmed" %in% colnames(fd_df) && "rt" %in% colnames(fd_df))
      fd_df$rtmed <- fd_df$rt
  }
  has_fd <- !is.null(fd_df)

  own_con <- is.null(con)
  if (own_con && !dir.exists(file_directory))
    cli::cli_abort("{.arg file_directory} does not exist: {.path {file_directory}}")

  # No sample produced a single pointer row anywhere: setup_connection()
  # would abort on the empty glob. Skip the DB connection entirely and
  # report every target as MissingMs1Signal instead - this degenerates to
  # the normal per-target sentinel-row path below rather than a hard
  # failure, so a target directory that legitimately has zero signal is a
  # labeled result, not a crash. A directory that doesn't exist at all
  # (above) is a usage error, not a data-quality one, and still aborts.
  if (own_con && !length(list.files(file_directory, pattern = "_pointers\\.parquet$"))) {
    return(data.table::rbindlist(lapply(targets, function(target) {
      mz_target <- if (has_fd && target %in% rownames(fd_df)) fd_df[target, "mzmed"] else NA_real_
      if (length(mz_target) != 1) mz_target <- NA_real_
      .deconv_error_row(target, "MissingMs1Signal", mz_target)
    }), fill = TRUE))
  }

  if (own_con) {
    con <- setup_connection(file_directory)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  if (report) {
    dir.create(dirname(report_path), showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(report_path, width = 11, height = 8)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  if (!is.null(mgf_dir)) dir.create(mgf_dir, showWarnings = FALSE, recursive = TRUE)

  total <- length(targets)
  fetch_chunks <- split(targets, ceiling(seq_along(targets) / fetch_chunk_size))

  pb_env <- environment()
  if (progress) {
    cli::cli_progress_bar(
      total  = total,
      format = "  {cli::pb_bar} {cli::pb_current}/{cli::pb_total}  {cli::pb_percent}  ETA {cli::pb_eta}  {cli::pb_status}",
      .envir = pb_env
    )
  }

  results <- data.table::rbindlist(lapply(fetch_chunks, function(chunk_ids) {
    raw <- batch_get_feature_data(chunk_ids, con, chunk_size = length(chunk_ids))
    if (nrow(raw) == 0) {
      if (progress) cli::cli_progress_update(inc = length(chunk_ids), .envir = pb_env)
      return(data.table::rbindlist(lapply(chunk_ids, function(target) {
        mz_target <- if (has_fd && target %in% rownames(fd_df)) fd_df[target, "mzmed"] else NA_real_
        if (length(mz_target) != 1) mz_target <- NA_real_
        .deconv_error_row(target, "MissingMs1Signal", mz_target)
      }), fill = TRUE))
    }
    by_feature <- split(raw, by = "feature")

    data.table::rbindlist(lapply(chunk_ids, function(target) {
      if (progress) cli::cli_progress_update(status = target, .envir = pb_env)

      dt <- by_feature[[target]]

      mz_target <- if (has_fd && target %in% rownames(fd_df)) fd_df[target, "mzmed"]
                   else if (!is.null(dt) && nrow(dt) > 0) dt$mz[1]
                   else NA_real_
      rt_target <- if (has_fd && target %in% rownames(fd_df)) fd_df[target, "rtmed"]
                   else NA_real_
      if (length(mz_target) != 1) mz_target <- NA_real_
      if (length(rt_target) != 1) rt_target <- NA_real_

      if (is.na(mz_target) || (require_rt && is.na(rt_target)))
        return(.deconv_error_row(target, "MissingTargetInfo", mz_target))

      if (is.null(dt) || nrow(dt) == 0)
        return(.deconv_error_row(target, "MissingMs1Signal", mz_target))

      stats <- tryCatch(
        process_pseudo_spec(dt, mz_target = mz_target,
                            grouper = grouper,
                            group_statistics_fn = group_statistics_fn,
                            buffer = buffer,
                            consensus_mz_fn = consensus_mz_fn),
        error = function(e) NULL
      )
      if (is.null(stats) || nrow(stats) == 0)
        return(.deconv_error_row(target, "NoAssociatedFragments", mz_target))

      pseudo_spec <- filter_fn(stats)
      if (is.null(pseudo_spec) || nrow(pseudo_spec) == 0)
        return(.deconv_error_row(target, "NoAssociatedFragments", mz_target))

      if (report) {
        gs <- attr(stats, "grouped_signal")
        plots <- Filter(Negate(is.null),
                        .pseudo_spec_plots(gs, stats, pseudo_spec, target))
        if (length(plots)) {
          if (requireNamespace("cowplot", quietly = TRUE)) {
            print(cowplot::plot_grid(plotlist = plots, ncol = 2))
          } else {
            for (p in plots) print(p)
          }
        }
      }

      if (!is.null(mgf_dir)) {
        write_simple_mgf(
          pseudo_spec, mz_target = mz_target, feature_id = target,
          path = file.path(mgf_dir, paste0(target, ".mgf")),
          rt = rt_target
        )
      }

      data.table::data.table(
        mz_g2           = pseudo_spec$fragment_mz,
        grouped_mz      = as.integer(pseudo_spec$group),
        mz              = pseudo_spec$fragment_mz,
        total_intensity = pseudo_spec$int,
        n_samples       = as.integer(pseudo_spec$n_samples),
        count           = as.integer(pseudo_spec$count),
        r               = pseudo_spec$correlation,
        t_stat          = pseudo_spec$correlation * sqrt(pseudo_spec$count - 2) /
         sqrt(1 - pseudo_spec$correlation^2),
        p_value         = pseudo_spec$p_value,
        fdr             = pseudo_spec$fdr,
        feature_id      = target,
        deconv_error    = NA_character_
      )
    }), fill = TRUE)
  }), fill = TRUE)

  if (progress) cli::cli_progress_done(.envir = pb_env)
  results
}
