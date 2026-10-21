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
#' Only `mz` is required; each grouper ignores the arguments it does not
#' use. [process_pseudo_spec()] calls the grouper **once per feature**,
#' on the MS2 centroids of all samples pooled together, passing
#' `mz` (fragment m/z), `intensity` (centroid intensity) and `target_mz`
#' (precursor m/z). `scan` and `ms1_int` are not passed there; they exist
#' for custom callers.
#'
#' Because samples are pooled, the tolerance has to cover between-sample
#' calibration drift as well as within-sample scatter.
#'
#' A grouper returns an integer vector, one group id per element of
#' `mz`, with no more distinct ids than distinct `mz` values. Check a
#' custom grouper with [validate_grouper()].
#'
#' @name grouper-contract
NULL


#' ppm-based UPGMA grouper (the package default)
#'
#' 1-D UPGMA on the sorted m/z values, cut so that no group is wider
#' (max m/z − min m/z) than `ppm` of the precursor m/z. Uses `mz` and
#' `target_mz`.
#'
#' @param ppm Maximum group width, in ppm of `target_mz`.
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
#' Agglomerative 1-D clustering where each centroid's m/z uncertainty
#' depends on its intensity: `sigma_i = mz_i * p0 * max(1, sqrt(Isat /
#' I_i)) / 1e6`. A group is kept only if every merge inside it is within
#' `k` sigma. Uses `mz` and `intensity`. See [fast_1D_iaac_grouping_cpp()].
#'
#' @param k Merge threshold in sigma. `k < 0` gives all singletons,
#'   `Inf` one group.
#' @param p0 Mass accuracy floor (ppm) reached at high intensity.
#' @param Isat Intensity at which the uncertainty reaches `p0`.
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
#' Exact optimal partition of the sorted m/z values: minimises the summed
#' per-group cost plus `beta` per group. Per-centroid uncertainty is the
#' same as in [iaac_grouper()]. Larger `beta` gives fewer, wider groups.
#' Uses `mz` and `intensity`. See [fast_1D_pelt_grouping_cpp()].
#'
#' @param beta Penalty per group; `<= 0` uses `2 * log(n)`.
#' @param p0 Mass accuracy floor (ppm) reached at high intensity.
#' @param Isat Intensity at which the uncertainty reaches `p0`.
#' @param max_seg Maximum centroids per group (`0` = no limit). Bounds
#'   run time at `O(n * max_seg)`; the result is optimal under the cap.
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
#' As [ppm_upgma()], but the maximum group width is a fixed value in Da.
#' Uses only `mz`.
#'
#' @param max_diff Maximum group width (max m/z − min m/z) in Da.
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
#' Group id is `floor(mz / bin_width) + 1`. Fast, but two close values
#' on either side of a bin edge end up in different groups. Uses only
#' `mz`.
#'
#' @param bin_width Bin width in Da.
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
#' Calls `grouper` on a small synthetic input and errors unless it
#' returns one id per input value and no more groups than distinct m/z
#' values. It does not check the type of the ids.
#'
#' @param grouper A function conforming to [grouper-contract].
#' @param mz_example Numeric vector to probe with.
#' @param intensity,target_mz,scan,ms1_int Probe values passed to
#'   `grouper`.
#' @return Invisibly `TRUE`; errors otherwise.
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
