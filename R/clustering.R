#' Grouper contract
#'
#' All groupers in this package share one signature so they are
#' drop-in replaceable in [clean_pseudo_spectra()] and
#' [process_pseudo_spec()]:
#'
#' \preformatted{
#'   function(mz, intensity = NULL, target_mz = NULL, scan = NULL, ms1_int = NULL)
#' }
#'
#' Only `mz` is required. Each grouper uses the subset of arguments it
#' needs and ignores the rest; callers pass whatever is available at
#' each grouping site (some columns are dropped by aggregation between
#' the two grouping passes in `process_pseudo_spec()`).
#'
#' A grouper must return an integer vector of the same length as `mz`
#' with no more unique values than `unique(mz)`; use [validate_grouper()]
#' to check a custom implementation.
#'
#' @name grouper-contract
#' @keywords internal
NULL


#' ppm-based UPGMA grouper (the package default)
#'
#' Returns a grouper closure that calls the package's C++ 1-D UPGMA
#' implementation with `maxDiff = 1e-6 * target_mz * ppm`.
#' Uses `mz` and `target_mz`; ignores `intensity`, `scan`, `ms1_int`.
#'
#' @param ppm Numeric. Tolerance in parts per million.
#' @return A grouper closure - see [grouper-contract].
#' @seealso [abs_upgma()], [iaac_grouper()], [pelt_grouper()], [validate_grouper()]
#' @export
#' @examples
#' g <- ppm_upgma(ppm = 30)
#' g(c(100.0, 100.001, 200.0), target_mz = 100)
ppm_upgma <- function(ppm = 30) {
  force(ppm)
  function(mz, intensity = NULL, target_mz = NULL, scan = NULL, ms1_int = NULL) {
    if (is.null(target_mz)) stop("ppm_upgma requires `target_mz`.")
    fast_1D_upgma_grouping_cpp(mz, maxDiff = 1e-6 * target_mz * ppm)
  }
}


#' Intensity-adaptive agglomerative clustering grouper
#'
#' Returns a grouper closure that calls the package's C++ 1-D IAAC
#' implementation. Uses `mz` and `intensity`; ignores `target_mz`,
#' `scan`, `ms1_int`.
#'
#' @param k Numeric. Number of clusters.
#' @param p0 Numeric. Initial parameter for clustering.
#' @param Isat Numeric. Saturation intensity.
#' @return A grouper closure - see [grouper-contract].
#' @seealso [ppm_upgma()], [abs_upgma()], [pelt_grouper()], [validate_grouper()]
#' @export
#' @examples
#' g <- iaac_grouper(k = 1, p0 = 8, Isat = 1e5)
#' g(c(100.0, 100.001, 200.0), intensity = c(1000, 800, 200))
iaac_grouper <- function(k = 1, p0 = 8, Isat = 1e5) {
  force(k); force(p0); force(Isat)
  function(mz, intensity = NULL, target_mz = NULL, scan = NULL, ms1_int = NULL) {
    if (is.null(intensity)) stop("iaac_grouper requires `intensity`.")
    fast_1D_iaac_grouping_cpp(mz, intensity = intensity,
                              k = k, p0 = p0, Isat = Isat)
  }
}


#' PELT (penalized optimal partitioning) grouper
#'
#' Returns a grouper closure that calls the package's C++ 1-D PELT
#' implementation. Unlike the greedy UPGMA / IAAC variants, PELT
#' globally minimises the total per-segment shot-noise cost plus a
#' `beta`-weighted segment count — the exact optimum over 1-D
#' partitions. Uses `mz` and `intensity`; ignores `target_mz`, `scan`,
#' `ms1_int`.
#'
#' Cost model matches [iaac_grouper()]: per-point sigma
#' `mz_i * p0 * max(1, sqrt(Isat / I_i)) / 1e6` (Da). Larger `beta`
#' means fewer, wider groups; `beta <= 0` selects the BIC-style default
#' `2 * log(n)`. `max_seg > 0` caps the largest allowed group size and
#' bounds the worst case at `O(n * max_seg)` while preserving optimality
#' among partitions obeying the cap.
#'
#' @param beta Numeric. Penalty per additional segment; `<= 0` picks
#'   `2 * log(n)`.
#' @param p0 Numeric. Floor mass accuracy in ppm.
#' @param Isat Numeric. Intensity at which per-point precision plateaus at `p0`.
#' @param max_seg Integer. Maximum group size; `0` = unlimited.
#' @return A grouper closure - see [grouper-contract].
#' @seealso [iaac_grouper()], [ppm_upgma()], [abs_upgma()], [validate_grouper()]
#' @export
#' @examples
#' g <- pelt_grouper(beta = -1, p0 = 5, Isat = 5e4)
#' g(c(100.000, 100.001, 100.500), intensity = c(1000, 800, 200))
pelt_grouper <- function(beta = -1.0, p0 = 5.0, Isat = 5e4, max_seg = 0L) {
  force(beta); force(p0); force(Isat); force(max_seg)
  function(mz, intensity = NULL, target_mz = NULL, scan = NULL, ms1_int = NULL) {
    if (is.null(intensity)) stop("pelt_grouper requires `intensity`.")
    fast_1D_pelt_grouping_cpp(mz, intensity = intensity,
                              beta = beta, p0 = p0, Isat = Isat,
                              max_seg = as.integer(max_seg))
  }
}


#' Absolute-tolerance UPGMA grouper
#'
#' Same backend as [ppm_upgma()] but with a fixed `max_diff` in Da.
#' Uses only `mz`.
#'
#' @param max_diff Numeric. Absolute m/z tolerance in Da.
#' @return A grouper closure - see [grouper-contract].
#' @export
#' @examples
#' g <- abs_upgma(max_diff = 0.005)
#' g(c(100.0, 100.001, 200.0))
abs_upgma <- function(max_diff = 0.005) {
  force(max_diff)
  function(mz, intensity = NULL, target_mz = NULL, scan = NULL, ms1_int = NULL) {
    fast_1D_upgma_grouping_cpp(mz, maxDiff = max_diff)
  }
}


#' Fixed-width binning grouper
#'
#' Assigns each m/z to a bin of width `bin_width` (Da) - fast,
#' deterministic. Bin-boundary values straddling a boundary get split
#' even when closer to each other than `bin_width`; for tolerant
#' grouping use [ppm_upgma()] / [abs_upgma()] instead. Uses only `mz`.
#'
#' @param bin_width Numeric. Bin width in Da (typical: 0.01, 0.05, 0.1).
#' @return A grouper closure - see [grouper-contract].
#' @export
#' @examples
#' g <- bin_grouper(bin_width = 0.01)
#' g(c(100.001, 100.009, 100.015))
bin_grouper <- function(bin_width = 0.01) {
  force(bin_width)
  function(mz, intensity = NULL, target_mz = NULL, scan = NULL, ms1_int = NULL) {
    as.integer(floor(mz / bin_width)) + 1L
  }
}


#' Sanity-check a custom grouper against the contract
#'
#' Probes a grouper with small example vectors and asserts the return
#' has one integer per input and no more groups than distinct m/z
#' values. Use once when wiring up a custom grouper; do not call inside
#' hot loops.
#'
#' @param grouper A function conforming to [grouper-contract].
#' @param mz_example Numeric vector to probe with.
#' @param intensity,target_mz,scan,ms1_int Probe values for groupers
#'   that require them; defaults are sensible synthetic inputs.
#' @return Invisibly `TRUE` when valid; throws otherwise.
#' @keywords internal
#' @export
validate_grouper <- function(grouper,
                             mz_example = c(100.000, 100.001, 100.500, 200.000),
                             intensity  = c(1000, 800, 200, 5000),
                             target_mz  = 100,
                             scan       = c(1L, 1L, 2L, 2L),
                             ms1_int    = c(5000, 5000, 4000, 4000)) {
  if (!is.function(grouper)) stop("`grouper` must be a function.")
  out <- grouper(
    mz        = mz_example,
    intensity = intensity,
    target_mz = target_mz,
    scan      = scan,
    ms1_int   = ms1_int
  )
  if (length(out) != length(mz_example)) {
    stop(sprintf("grouper returned length %d, expected %d.",
                 length(out), length(mz_example)))
  }
  if (length(unique(out)) > length(unique(mz_example))) {
    stop("grouper produced more unique groups than input m/z values.")
  }
  invisible(TRUE)
}
