#' Statistics contract
#'
#' All group-statistics closures in this package share one signature so they
#' are drop-in replaceable in [process_pseudo_spec()] and
#' [clean_pseudo_spectra()]:
#'
#' \preformatted{
#'   function(dt) -> per-group summary data.table
#' }
#'
#' `dt` is the raw per-feature data.table produced by
#' [batch_get_feature_data()] (columns `sample, feature, mz, fragment_mz,
#' int, ms1_scan, ms2_scan, ms1_int`) after the apex window has been
#' applied, the active grouper has added a `group` column, and
#' [process_pseudo_spec()] has added `group_mz` — the consensus m/z, one
#' constant value per group.
#'
#' A statistics closure must **carry `group_mz` through as `fragment_mz`**
#' rather than deriving m/z itself; m/z is settled before statistics run.
#' It aggregates and returns **one row per group** with the canonical
#' schema
#'
#' \preformatted{
#'   group, fragment_mz, mz, int, ms1_int,
#'   correlation, p_value, fdr, count, n_samples
#' }
#'
#' plus any extra per-group columns it needs (e.g. `apex_shift_med`,
#' `apex_shift_iqr` from [apex_correlation_statistics()]). No thresholding
#' is applied — that is the job of a [filter-contract] closure downstream.
#'
#' It should attach the per-`(group, sample, ms1_scan)` aggregation as an
#' attribute named `"grouped_signal"` so [process_pseudo_spec()] and the
#' diagnostic plots can reuse it without recomputing.
#'
#' Return an empty data.table (same columns) when nothing survives
#' `min_count`.
#'
#' @name stats-contract
#' @keywords internal
NULL


#' Filter contract
#'
#' All filter closures in this package share one signature so they are
#' drop-in replaceable in [clean_pseudo_spectra()]:
#'
#' \preformatted{
#'   function(stats) -> row-subset of stats
#' }
#'
#' `stats` is the per-group summary produced by a [stats-contract] closure.
#' A filter returns a row-subset of that summary (same columns) — no
#' aggregation, no new statistics.
#'
#' Return an empty data.table (same columns) when nothing survives.
#'
#' @name filter-contract
#' @keywords internal
NULL


#' Aggregate + Pearson correlation + BH FDR (shared statistics internals)
#'
#' Aggregates `dt` to `(group, sample, ms1_scan)`, prunes groups below
#' `min_count`, computes per-group Pearson r and p-value of summed MS2
#' intensity against MS1 intensity, and BH-adjusts p-values across
#' surviving groups. Returns the standard per-group summary
#' (unfiltered — thresholds are the filter closure's job). The
#' `(group, sample, ms1_scan)` aggregation is attached as attribute
#' `"grouped_signal"` so downstream code (extra-statistic add-ons,
#' diagnostic plots) can reuse it without recomputing.
#'
#' m/z is not computed here. [process_pseudo_spec()] has already written
#' one consensus value per group into `group_mz`, and this function only
#' carries it through as `fragment_mz`. `dt` arriving without that column
#' (a direct notebook call on [get_raw_spec()] output) gets it filled in
#' first, so the estimate is the same either way. The `fragment_mz` inside
#' `grouped_signal` is a plain per-scan-cell mean kept for diagnostics —
#' nothing reads it for the estimate.
#'
#' @keywords internal
#' @noRd
.correlations_from_dt <- function(dt, min_count) {
  if (!"group_mz" %in% names(dt))
    dt[, group_mz := .consensus_mz(fragment_mz, int), by = group]

  grouped_signal <- dt[, .(
    int         = sum(int),
    ms1_int     = mean(ms1_int),
    fragment_mz = mean(fragment_mz),
    mz          = mean(mz),
    group_mz    = group_mz[1]
  ), by = .(group, sample, ms1_scan)]

  group_counts <- grouped_signal[, .(count = .N), by = group][count > min_count]
  if (nrow(group_counts) == 0) {
    out <- .empty_correlations()
    data.table::setattr(out, "grouped_signal", grouped_signal[0])
    return(out)
  }

  cor_results <- grouped_signal[group %in% group_counts$group, {
    n <- .N
    if (n < 3L) {
      list(p_value = NA_real_, correlation = NA_real_)
    } else {
      xc <- int - mean(int)
      yc <- ms1_int - mean(ms1_int)
      denom <- sqrt(sum(xc * xc) * sum(yc * yc))
      if (denom == 0) {
        list(p_value = NA_real_, correlation = NA_real_)
      } else {
        r <- max(-1, min(1, sum(xc * yc) / denom))
        if (abs(r) >= 1) {
          list(p_value = 0, correlation = r)
        } else {
          tval <- r * sqrt((n - 2L) / (1 - r * r))
          list(p_value = 2 * stats::pt(-abs(tval), df = n - 2L), correlation = r)
        }
      }
    }
  }, by = group]

  group_summary <- grouped_signal[group %in% group_counts$group, .(
    fragment_mz = group_mz[1],
    mz          = mean(mz),
    int         = mean(int),
    ms1_int     = mean(ms1_int),
    n_samples   = data.table::uniqueN(sample)
  ), by = group]

  out <- merge(group_summary, cor_results, by = "group")
  out <- merge(out, group_counts, by = "group")
  out[, fdr := stats::p.adjust(p_value, method = "fdr")]
  data.table::setattr(out, "grouped_signal", grouped_signal)
  out
}

#' Empty per-group summary with the canonical schema.
#' @keywords internal
#' @noRd
.empty_correlations <- function() {
  data.table::data.table(
    group       = integer(0),
    fragment_mz = numeric(0),
    mz          = numeric(0),
    int         = numeric(0),
    ms1_int     = numeric(0),
    n_samples   = integer(0),
    correlation = numeric(0),
    p_value     = numeric(0),
    count       = integer(0),
    fdr         = numeric(0)
  )
}


#' Per-group correlation statistics (the package default)
#'
#' Returns a statistics closure that aggregates `dt` to
#' `(group, sample, ms1_scan)`, drops groups with fewer than `min_count`
#' observations, and computes per-group Pearson correlation between
#' summed MS2 intensity and MS1 intensity, its p-value, and the
#' BH-adjusted FDR. No thresholding is applied — see
#' [correlation_filtering()] for the paired filter.
#'
#' The reported `fragment_mz` is the per-group consensus m/z computed
#' upstream by [process_pseudo_spec()]; this closure carries it through
#' and does not recompute it.
#'
#' Conforms to [stats-contract].
#'
#' @param min_count Minimum `(sample, ms1_scan)` observation count per
#'   group before Pearson r is computed.
#' @return A statistics closure — see [stats-contract].
#' @seealso [apex_correlation_statistics()], [correlation_filtering()]
#' @export
correlation_statistics <- function(min_count = 10L) {
  force(min_count)
  function(dt) {
    if (nrow(dt) == 0) return(.empty_correlations())
    .correlations_from_dt(dt, min_count = min_count)
  }
}


#' Per-group correlation + apex-agreement statistics
#'
#' As [correlation_statistics()], plus two extra per-group columns
#' summarising the intra-sample MS1 <-> MS2 apex-centroid shift across
#' samples:
#'
#' \itemize{
#'   \item `apex_shift_med` — median across samples of the intensity-weighted
#'     `(m2 - m1)` scan-index centroid difference.
#'   \item `apex_shift_iqr` — IQR of the same across samples.
#' }
#'
#' MS1 centroids are computed on the deduplicated per-`(sample, ms1_scan)`
#' MS1 intensity trace so each scan is weighted once by `ms1_int`, not
#' `ms1_int x N_fragments_at_that_scan` (the earlier all-rows form biased
#' the centroid toward groups whose fragments cluster at the apex).
#'
#' Conforms to [stats-contract]. Pair with [apex_correlation_filtering()],
#' which additionally thresholds on the two apex-shift columns.
#'
#' @inheritParams correlation_statistics
#' @return A statistics closure — see [stats-contract].
#' @seealso [correlation_statistics()], [apex_correlation_filtering()]
#' @export
apex_correlation_statistics <- function(min_count = 10L) {
  force(min_count)
  function(dt) {
    if (nrow(dt) == 0) return(.empty_correlations())
    out <- .correlations_from_dt(dt, min_count = min_count)
    if (nrow(out) == 0) {
      out[, `:=`(apex_shift_med = numeric(0),
                 apex_shift_iqr = numeric(0))]
      return(out)
    }

    ms1_by_sample <- unique(dt[, .(sample, ms1_scan, ms1_int)])[
      , .(m1 = sum(ms1_scan * ms1_int) / sum(ms1_int)), by = sample]
    ms2_by_group <- dt[!is.na(fragment_mz),
      .(m2 = sum(ms1_scan * int) / sum(int)), by = .(sample, group)]
    apex_per_sample <- merge(ms2_by_group, ms1_by_sample, by = "sample")
    apex_per_sample[, apex_shift := m2 - m1]

    apex_per_group <- apex_per_sample[, .(
      apex_shift_med = stats::median(apex_shift),
      apex_shift_iqr = stats::IQR(apex_shift)
    ), by = group]

    gs <- attr(out, "grouped_signal")
    out <- merge(out, apex_per_group, by = "group", all.x = TRUE)
    data.table::setattr(out, "grouped_signal", gs)
    out
  }
}


#' Correlation + BH FDR filter (the package default)
#'
#' Returns a filter closure that keeps rows whose Pearson correlation is
#' at least `cor_thresh` and whose BH-adjusted FDR is at most
#' `fdr_thresh`. Operates on the per-group summary produced by a
#' [stats-contract] closure (typically [correlation_statistics()]).
#'
#' Conforms to [filter-contract].
#'
#' @param cor_thresh Minimum Pearson correlation to keep.
#' @param fdr_thresh Maximum BH-adjusted FDR to keep.
#' @return A filter closure — see [filter-contract].
#' @seealso [apex_correlation_filtering()], [correlation_statistics()]
#' @export
correlation_filtering <- function(cor_thresh = 0.3, fdr_thresh = 0.01) {
  force(cor_thresh); force(fdr_thresh)
  function(stats) {
    if (nrow(stats) == 0) return(stats)
    stats[correlation >= cor_thresh & fdr <= fdr_thresh]
  }
}


#' Correlation + FDR + apex-agreement filter
#'
#' As [correlation_filtering()], additionally requiring
#' `abs(apex_shift_med) <= max_apex_shift_med` and
#' `apex_shift_iqr <= max_apex_shift_iqr`. Requires the two apex-shift
#' columns on its input — pair with [apex_correlation_statistics()]
#' upstream; a helpful error is raised if those columns are missing.
#'
#' Conforms to [filter-contract].
#'
#' @inheritParams correlation_filtering
#' @param max_apex_shift_med Reject if `abs(median(apex_shift))` across
#'   samples exceeds this (in scans).
#' @param max_apex_shift_iqr Reject if `IQR(apex_shift)` across samples
#'   exceeds this (in scans).
#' @return A filter closure — see [filter-contract].
#' @seealso [correlation_filtering()], [apex_correlation_statistics()]
#' @export
apex_correlation_filtering <- function(cor_thresh         = 0.3,
                                       fdr_thresh         = 0.01,
                                       max_apex_shift_med = 0.5,
                                       max_apex_shift_iqr = 0.75) {
  force(cor_thresh); force(fdr_thresh)
  force(max_apex_shift_med); force(max_apex_shift_iqr)
  function(stats) {
    if (nrow(stats) == 0) return(stats)
    required <- c("apex_shift_med", "apex_shift_iqr")
    missing_cols <- setdiff(required, colnames(stats))
    if (length(missing_cols)) {
      stop("apex_correlation_filtering requires column(s) [",
           paste(missing_cols, collapse = ", "),
           "]; use apex_correlation_statistics() upstream.")
    }
    stats[
      correlation           >= cor_thresh &
      fdr                   <= fdr_thresh &
      abs(apex_shift_med)   <= max_apex_shift_med &
      apex_shift_iqr        <= max_apex_shift_iqr
    ]
  }
}
