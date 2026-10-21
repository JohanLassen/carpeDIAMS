#' Open a DuckDB connection with views over pool/pointer parquets
#'
#' Opens an in-memory DuckDB with views `pointers` and `pool` over all
#' `*_pointers.parquet` / `*_pool.parquet` files in `file_directory`.
#' DuckDB spills to `<file_directory>/.duckdb_tmp`, with a fixed
#' `memory_limit = '128GB'` and `threads = 2`. Errors if there are no
#' pointer files. Close it with `DBI::dbDisconnect(con, shutdown = TRUE)`;
#' [clean_pseudo_spectra()] does not close a connection you pass in.
#'
#' @param file_directory Directory with the per-sample parquets.
#' @return A DuckDB connection.
#' @keywords internal
#' @export
setup_connection <- function(file_directory) {
  # Resolve to absolute path — DuckDB may not share R's working directory
  abs_dir <- normalizePath(file_directory, mustWork = TRUE)
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
    "CREATE OR REPLACE VIEW pointers AS SELECT * FROM read_parquet('%s/*_pointers.parquet', union_by_name = true)",
    abs_dir
  ))
  DBI::dbExecute(con, sprintf(
    "CREATE OR REPLACE VIEW pool AS SELECT * FROM read_parquet('%s/*_pool.parquet', union_by_name = true)",
    abs_dir
  ))
  DBI::dbExecute(con, sprintf("SET temp_directory = '%s'", tmp_dir))
  DBI::dbExecute(con, "SET memory_limit = '128GB'")   # below the SLURM allocation
  DBI::dbExecute(con, "SET threads = 2")
  con
}

#' Fetch raw MS2 + MS1 rows for a vector of features from DuckDB
#'
#' Pairs each MS1 pointer row with the MS2 centroids of the **next scan**
#' (`pool.scan = pointers.scan + 1`, same sample), keeping fragments with
#' m/z at most precursor + 2. MS1 rows with no matching MS2 centroid are
#' kept with `NA` fragment columns.
#'
#' @param feature_ids Character vector of feature ids.
#' @param con A connection from [setup_connection()].
#' @param chunk_size Features per query.
#' @return `data.table`, one row per (MS1 scan, MS2 centroid):
#'   `sample, feature, mz` (precursor), `fragment_mz, int, ms1_scan,
#'   ms2_scan, ms1_int, apex_scan, ms1_rt, ms2_rt`. `apex_scan` is the
#'   apex found upstream (xcms peak or nearest scan to the target RT);
#'   the RT columns are `NULL`/`NA` for parquets written without `rt`.
#' @keywords internal
#' @export
batch_get_feature_data <- function(feature_ids, con, chunk_size = 500) {
  has_rt <- "rt" %in% DBI::dbListFields(con, "pointers") && "rt" %in% DBI::dbListFields(con, "pool")
  rt_select <- if (has_rt) "p.rt AS ms1_rt,\n        m.rt AS ms2_rt" else "NULL AS ms1_rt,\n        NULL AS ms2_rt"

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
        p.ms1_int,
        p.apex_scan,
        %s
      FROM feat_pointers p
      LEFT JOIN pool m ON (
        p.scan + 1 = m.scan AND
        p.sample = m.sample
      )
      WHERE (m.mz <= p.mz_target + 2 OR m.mz IS NULL)
      ORDER BY p.feature, p.sample, p.scan
    ", id_list, rt_select)
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

#' Per-sample FWHM crossing scans, apex-anchored (not window-max-anchored)
#'
#' The half-max threshold is fixed at the intensity *at the apex scan*,
#' never at the trace's own max - a co-eluting neighbor peak elsewhere in
#' the (already apex-windowed-by-nothing-yet, i.e. still-full) trace can
#' have higher intensity than the true feature's own apex, and using the
#' trace's max as the reference would silently pick the neighbor's
#' half-height instead. Walks outward from the apex and stops at the
#' FIRST sub-threshold scan on each side (a literal FWHM, not a smoothed
#' or interpolated one); if no crossing is found before the edge of
#' whatever data was actually extracted upstream, clamps to that edge
#' rather than erroring, extrapolating, or falling back to any other
#' width - deliberately asymmetric in how forgiving it is between the two
#' data routes (see `08_carpeDIAMS/carpeDIAMS_optimization`'s Stage E
#' notes: the xcms route's extraction window has generous peak-shape-aware
#' slack, so a clamp is rare; the target-list route's fixed rt_width does
#' not, so a clamp there more often means "hit the extraction boundary",
#' not "found the true half-max").
#'
#' `scan`/`int` may contain the same `(scan, int)` pair repeated many
#' times (one row per MS2 fragment sharing that MS1 scan) - deduplicated
#' before the walk, since duplicates would not change the FWHM crossing
#' point but would make the row-index walk meaningless.
#'
#' @param scan Numeric/integer MS1 scan numbers (may contain duplicates).
#' @param int MS1 intensity at each `scan` (same length as `scan`).
#' @param apex_scan Numeric. The reference apex scan; not required to be
#'   an exact match in `scan` (falls back to the nearest one present).
#' @return `list(lo = <scan>, hi = <scan>)`.
#' @keywords internal
#' @noRd
.fwhm_bounds <- function(scan, int, apex_scan) {
  tr <- unique(data.table::data.table(scan = scan, int = int))
  data.table::setorder(tr, scan)
  apex_row <- which.min(abs(tr$scan - apex_scan))
  half <- tr$int[apex_row] / 2

  lo <- tr$scan[apex_row]
  if (apex_row > 1L) {
    for (i in (apex_row - 1L):1L) {
      if (tr$int[i] < half) break
      lo <- tr$scan[i]
    }
  }
  hi <- tr$scan[apex_row]
  if (apex_row < nrow(tr)) {
    for (i in (apex_row + 1L):nrow(tr)) {
      if (tr$int[i] < half) break
      hi <- tr$scan[i]
    }
  }
  list(lo = lo, hi = hi)
}

#' Apply `.fwhm_bounds()` per sample and filter `dt` to the resulting window
#'
#' Note: `dt` at this point has already had fragment-less pointer rows
#' dropped (`process_pseudo_spec()`'s `!is.na(int) & !is.na(fragment_mz)`
#' filter, applied before this is ever called) - so the MS1 trace this
#' walks is "scans that had at least one detected fragment", not the
#' complete raw MS1 trace. This matches `apex_source = "recompute"`'s
#' existing (pre-existing, unchanged) behavior exactly, so defaults stay
#' bit-identical, but it does mean a genuinely fragment-free scan near a
#' true half-max crossing is invisible to the walk - a real, inherited
#' limitation of computing the window downstream of fragment matching,
#' not something specific to FWHM. Left as-is deliberately: fixing it
#' would change `apex_source = "recompute"`'s pinned regression numbers
#' too, which is a separate, bigger decision than adding FWHM windowing.
#' @keywords internal
#' @noRd
.fwhm_window <- function(dt) {
  bounds <- dt[, .fwhm_bounds(ms1_scan, ms1_int, apex[1L]), by = sample]
  dt <- dt[bounds, on = "sample"]
  dt <- dt[ms1_scan >= lo & ms1_scan <= hi]
  dt[, c("lo", "hi") := NULL]
  dt
}

#' Build a feature's per-group statistics table from raw MS1+MS2 data.table
#'
#' Turns one feature's raw rows into a per-group statistics table:
#'
#' 1. Drop rows without a fragment.
#' 2. Find each sample's apex scan (`apex_source`) and keep the scans
#'    around it (`window_mode`).
#' 3. Group the fragment m/z of all samples together with `grouper`.
#' 4. Compute one consensus m/z per group (`consensus_mz_fn`).
#' 5. Score the groups with `group_statistics_fn`.
#'
#' No filtering is applied; [clean_pseudo_spectra()] adds that.
#'
#' @param dt One feature's rows from [batch_get_feature_data()].
#' @param mz_target Precursor m/z, passed to the grouper as `target_mz`.
#' @param grouper A grouper closure ([grouper-contract]).
#' @param group_statistics_fn A statistics closure ([stats-contract]).
#' @param buffer With `window_mode = "buffer"`, keep scans with
#'   `abs(ms1_scan - apex) < buffer`.
#' @param consensus_mz_fn A consensus-mz closure
#'   ([consensus-mz-contract]); `NULL` uses [consensus_mz_weighted()].
#' @param apex_source `"recompute"`: the scan with the highest MS1
#'   intensity in `dt`. `"upstream"`: `dt$apex_scan`; samples where it is
#'   `NA` are dropped with a warning.
#' @param window_mode `"buffer"`: fixed scan window (see `buffer`).
#'   `"fwhm"`: walk out from the apex until MS1 intensity drops below half
#'   the apex intensity, within the scans present in `dt`.
#' @return The statistics table; zero rows if nothing is left at any step.
#' @export
process_pseudo_spec <- function(
  dt,
  mz_target,
  grouper             = ppm_upgma(ppm = 30),
  group_statistics_fn = correlation_statistics(),
  buffer              = 50L,
  consensus_mz_fn     = NULL,
  apex_source         = c("recompute", "upstream"),
  window_mode         = c("buffer", "fwhm")
) {
  apex_source <- match.arg(apex_source)
  window_mode <- match.arg(window_mode)
  if (is.null(consensus_mz_fn)) consensus_mz_fn <- consensus_mz_weighted()
  if (nrow(dt) == 0) return(.empty_correlations())
  dt <- dt[!is.na(int) & !is.na(fragment_mz)]
  if (nrow(dt) == 0) return(.empty_correlations())

  if (apex_source == "recompute") {
    dt[, apex := ms1_scan[which.max(ms1_int)], by = sample]
  } else {
    if (!"apex_scan" %in% names(dt))
      cli::cli_abort(c(
        "{.arg apex_source} = \"upstream\" requires an {.field apex_scan} column in {.arg dt}.",
        "i" = "Select {.code p.apex_scan} in batch_get_feature_data()'s SQL, or pass it through by hand."
      ))
    na_samples <- dt[, all(is.na(apex_scan)), by = sample][V1 == TRUE, sample]
    if (length(na_samples))
      cli::cli_warn("Dropping sample{?s} {.val {na_samples}} - apex_scan is NA (didn't survive the MS2 join).")
    dt <- dt[!(sample %in% na_samples)]
    if (nrow(dt) == 0) return(.empty_correlations())
    dt[, apex := apex_scan[1L], by = sample]
  }

  if (window_mode == "buffer") {
    dt <- dt[abs(ms1_scan - apex) < buffer]
  } else {
    dt <- .fwhm_window(dt)
  }
  if (nrow(dt) == 0) return(.empty_correlations())

  dt[, group := grouper(mz        = fragment_mz,
                        intensity = int,
                        target_mz = mz_target)]
  dt[, group_mz := consensus_mz_fn(fragment_mz, int), by = group]

  group_statistics_fn(dt)
}

#' Legacy per-group thresholder (superseded by [correlation_filtering()])
#'
#' Not exported. Use `correlation_filtering(cor_thresh, fdr_thresh)`.
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
#' Runs [process_pseudo_spec()] and `filter_fn` for every target, reading
#' the per-sample parquets in `file_directory` through DuckDB. This is
#' the work behind the `filter` and `finalize` stages.
#'
#' Every target gives either fragment rows or exactly one error row with
#' `deconv_error` set to:
#'
#' * `"MissingTargetInfo"`: no precursor m/z (or no RT when
#'   `require_rt = TRUE`).
#' * `"MissingMs1Signal"`: no pointer rows for the target, or no pointer
#'   parquets at all.
#' * `"NoAssociatedFragments"`: no group survived statistics or filtering,
#'   or [process_pseudo_spec()] raised an error (errors are not shown).
#'
#' @param targets Character vector of feature ids.
#' @param file_directory Directory with the per-sample parquets.
#' @param fd Optional table with `mzmed` and `rtmed` (or `mz`/`rt`), with
#'   rownames or a `feature_id` column matching `targets`. Without it the
#'   precursor m/z comes from the parquets and RT is unknown.
#' @param con Optional open connection from [setup_connection()]; it is
#'   left open.
#' @param grouper A grouper closure ([grouper-contract]).
#' @param buffer Passed to [process_pseudo_spec()].
#' @param fdr_thresh,cor_thresh Used only to build the default `filter_fn`.
#' @param min_count Used only to build the default `group_statistics_fn`.
#' @param group_statistics_fn A statistics closure ([stats-contract]);
#'   `NULL` uses `correlation_statistics(min_count)`.
#' @param filter_fn A filter closure ([filter-contract]); `NULL` uses
#'   `correlation_filtering(cor_thresh, fdr_thresh)`.
#' @param consensus_mz_fn,apex_source,window_mode Passed to
#'   [process_pseudo_spec()].
#' @param require_rt If `TRUE`, a target without RT in `fd` is a
#'   `"MissingTargetInfo"` error. The target-list route sets this.
#' @param batch_idx,n_batches Process only slice `batch_idx` (1-based)
#'   of `targets` split into `n_batches` contiguous chunks. A slice past
#'   the end returns an empty table.
#' @param report,report_path Write diagnostic plots per feature to a PDF
#'   (needs ggplot2). **On by default**; the `process_*()` functions turn
#'   it off.
#' @param mgf_dir If set, also write one MGF file per feature here.
#' @param fetch_chunk_size Features per DuckDB query.
#' @param progress Show a progress bar.
#' @return `data.table` with one row per kept fragment: `mz_g2` and `mz`
#'   (both the consensus fragment m/z), `grouped_mz` (group id),
#'   `total_intensity` (mean per-cell intensity, not a sum), `n_samples`,
#'   `count`, `r`, `t_stat`, `p_value`, `fdr`, `feature_id`,
#'   `deconv_error` (`NA` on real rows). Column names are what
#'   [export_mgf()] reads.
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
  apex_source         = c("recompute", "upstream"),
  window_mode         = c("buffer", "fwhm"),
  require_rt          = FALSE,
  batch_idx           = NULL,
  n_batches           = 1L,
  report              = TRUE,
  report_path         = "pseudo_specs_report.pdf",
  mgf_dir             = NULL,
  fetch_chunk_size    = 50L,
  progress            = TRUE
) {
  apex_source <- match.arg(apex_source)
  window_mode <- match.arg(window_mode)
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
                            consensus_mz_fn = consensus_mz_fn,
                            apex_source = apex_source,
                            window_mode = window_mode),
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
