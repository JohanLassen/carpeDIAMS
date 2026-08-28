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



#' Build a feature's per-group statistics table from raw MS1+MS2 data.table
#'
#' Per-sample apex window (buffer) → first grouper pass on raw
#' `fragment_mz` → aggregation to per-first-pass-group mean m/z → second
#' grouper pass (merges neighbours that scan-to-scan m/z jitter would
#' otherwise split) → statistics closure returns the full per-group
#' summary. No thresholding — see [clean_pseudo_spectra()] for the
#' composition with a [filter-contract] closure.
#'
#' The `(group, sample, ms1_scan)` aggregation computed inside the
#' statistics closure is attached to the return value as attribute
#' `"grouped_signal"` so downstream code (diagnostic plots) can reuse it
#' without a second groupby.
#'
#' @param dt Output of [batch_get_feature_data()] restricted to one feature.
#' @param mz_target Numeric. Precursor m/z (drives the ppm tolerance).
#' @param grouper A grouper closure - see [grouper-contract].
#' @param group_statistics_fn A statistics closure - see [stats-contract].
#' @param buffer Integer. Scan window around each sample's MS1 apex.
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
  buffer              = 50L
) {
  if (nrow(dt) == 0) return(.empty_correlations())
  dt <- dt[!is.na(int) & !is.na(fragment_mz)]
  if (nrow(dt) == 0) return(.empty_correlations())

  dt[, apex := ms1_scan[which.max(ms1_int)], by = sample]
  dt <- dt[abs(ms1_scan - apex) < buffer]
  if (nrow(dt) == 0) return(.empty_correlations())

  # first grouping pass — on raw per-scan fragment_mz
  dt[, group := grouper(mz        = fragment_mz,
                        intensity = int,
                        target_mz = mz_target,
                        scan      = ms2_scan,
                        ms1_int   = ms1_int)]

  # second pass — on per-group mean m/z; neighbours split by
  # scan-to-scan m/z jitter get re-merged. The aggregation drops
  # scan/ms1_int, so groupers that require those are single-pass by
  # necessity — the second call still receives them as NULL.
  group_means <- dt[, .(fragment_mz = mean(fragment_mz),
                        intensity   = sum(int)),
                    by = .(group_first = group)]
  group_means[, group := grouper(mz        = fragment_mz,
                                 intensity = intensity,
                                 target_mz = mz_target)]
  dt[group_means, on = .(group = group_first), group := i.group]

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

#' Clean pseudo-spectra for many features (main batch entry)
#'
#' Orchestrates: connect → fetch → per-feature [process_pseudo_spec()]
#' (grouper + `group_statistics_fn`) → `filter_fn` → optional plot / MGF.
#' Returns a `data.table` in the schema downstream `03_make_mgf.R` expects.
#'
#' @param targets Character vector of feature ids.
#' @param file_directory Path to pool/pointer parquets.
#' @param fd Optional `featureDefinitions()`-shaped frame; supplies
#'   `mzmed` and `rtmed` per target.
#' @param con Optional pre-opened connection (else opened + closed internally).
#' @param buffer,fdr_thresh,cor_thresh,min_count Tunables passed downstream.
#'   `min_count` is consumed by the default `group_statistics_fn`;
#'   `cor_thresh` / `fdr_thresh` are consumed by the default `filter_fn`.
#'   A custom `group_statistics_fn` / `filter_fn` may ignore them.
#' @param grouper A grouper closure - see [grouper-contract]. Called
#'   twice inside [process_pseudo_spec()] (first on raw fragment m/z,
#'   then on per-group means to re-merge neighbours split by
#'   scan-to-scan m/z jitter).
#' @param group_statistics_fn A statistics closure conforming to
#'   [stats-contract] — receives the raw per-feature dt (with `group`
#'   column) and returns the per-group summary (no thresholding).
#'   `NULL` (default) uses [correlation_statistics()] with `min_count`.
#'   Swap in [apex_correlation_statistics()] to add apex-shift columns.
#' @param filter_fn A filter closure conforming to [filter-contract] —
#'   receives the per-group stats table and returns a row-subset.
#'   `NULL` (default) uses [correlation_filtering()] with
#'   `cor_thresh`/`fdr_thresh`. Swap in [apex_correlation_filtering()]
#'   for apex-aware thresholding (pair with
#'   [apex_correlation_statistics()]).
#' @param batch_idx,n_batches Optional slice: process only batch
#'   `batch_idx` of `n_batches` equal contiguous chunks of `targets`.
#'   `batch_idx = NULL` (default) processes every target in-process.
#'   Set both for the DAG-batched `stage = "filter"` path.
#' @param report,report_path PDF diagnostic report on/off + path.
#' @param mgf_dir Optional directory for per-feature MGFs.
#' @param fetch_chunk_size Features per DuckDB query.
#' @param progress Print a per-feature progress counter.
#' @return A `data.table` matching the legacy 02_filter_features schema.
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
  own_con <- is.null(con)
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
      return(NULL)
    }
    by_feature <- split(raw, by = "feature")

    data.table::rbindlist(lapply(chunk_ids, function(target) {
      if (progress) cli::cli_progress_update(status = target, .envir = pb_env)

      dt <- by_feature[[target]]
      if (is.null(dt) || nrow(dt) == 0) return(NULL)

      mz_target <- if (has_fd && target %in% rownames(fd_df)) fd_df[target, "mzmed"]
                   else dt$mz[1]
      if (length(mz_target) != 1 || is.na(mz_target)) return(NULL)
      rt_target <- if (has_fd && target %in% rownames(fd_df)) fd_df[target, "rtmed"]
                   else NA_real_

      stats <- tryCatch(
        process_pseudo_spec(dt, mz_target = mz_target,
                            grouper = grouper,
                            group_statistics_fn = group_statistics_fn,
                            buffer = buffer),
        error = function(e) NULL
      )
      if (is.null(stats) || nrow(stats) == 0) return(NULL)

      pseudo_spec <- filter_fn(stats)
      if (is.null(pseudo_spec) || nrow(pseudo_spec) == 0) return(NULL)

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
        feature_id      = target
      )
    }), fill = TRUE)
  }), fill = TRUE)

  if (progress) cli::cli_progress_done(.envir = pb_env)
  results
}
