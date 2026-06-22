#' ppm-based UPGMA grouper (the package default)
#'
#' Returns a grouper closure that calls the package's C++ 1-D UPGMA
#' implementation with `maxDiff = 1e-6 * target_mz * ppm`.
#'
#' @param ppm Numeric. Tolerance in parts per million.
#' @return A `function(mz, target_mz)` satisfying the grouper contract.
#' @seealso [abs_upgma()], [validate_grouper()], `fast_1D_upgma_grouping_cpp()`
#' @export
#' @examples
#' g <- ppm_upgma(ppm = 30)
#' g(c(100.0, 100.001, 200.0), target_mz = 100)
ppm_upgma <- function(ppm = 30) {
force(ppm)
function(mz, target_mz) {
    fast_1D_upgma_grouping_cpp(mz, maxDiff = 1e-6 * target_mz * ppm)
}
}

#' ppm-based UPGMA grouper (the package default)
#'
#' Returns a grouper closure that calls the package's C++ 1-D UPGMA
#' implementation with `maxDiff = 1e-6 * target_mz * ppm`.
#'
#' @param ppm Numeric. Tolerance in parts per million.
#' @return A `function(mz, target_mz)` satisfying the grouper contract.
#' @seealso [abs_upgma()], [validate_grouper()], `fast_1D_upgma_grouping_cpp()`
#' @export
#' @examples
#' g <- ppm_upgma(ppm = 30)
#' g(c(100.0, 100.001, 200.0), target_mz = 100)
ppm_upgma <- function(ppm = 30) {
force(ppm)
function(mz, target_mz) {
    fast_1D_upgma_grouping_cpp(mz, maxDiff = 1e-6 * target_mz * ppm)
}
}

#' Intensity adaptive agglomerative clustering grouper 
#'
#' Returns a grouper closure that calls the package's C++ 1-D IAAC
#' implementation with `k`, `p0`, and `Isat` parameters.
#'
#' @param k Numeric. Number of clusters.
#' @param p0 Numeric. Initial parameter for clustering.
#' @param Isat Numeric. Saturation intensity.
#' @return A `function(mz)` satisfying the grouper contract.
#' @seealso [abs_upgma()], [validate_grouper()], `fast_1D_iaac_grouping_cpp()`
#' @export
#' @examples
#' g <- iaac_grouper(k = 1, p0 = 8, Isat = 1e5)
#' g(c(100.0, 100.001, 200.0))
iaac_grouper <- function(k = 1, p0 = 8, Isat = 1e5) {
force(k)
force(p0)
force(Isat)
function(mz) {
    fast_1D_iaac_grouping_cpp(mz, k = k, p0 = p0, Isat = Isat)
}
}



#' Absolute-tolerance UPGMA grouper
#'
#' Same backend as [ppm_upgma()], but uses a fixed `max_diff` in Da and
#' ignores `target_mz`.
#'
#' @param max_diff Numeric. Absolute m/z tolerance in Da.
#' @return A `function(mz, target_mz)` satisfying the grouper contract.
#' @export
#' @examples
#' g <- abs_upgma(max_diff = 0.005)
#' g(c(100.0, 100.001, 200.0), target_mz = NA)
abs_upgma <- function(max_diff = 0.005) {
force(max_diff)
function(mz, target_mz = NULL) {
    fast_1D_upgma_grouping_cpp(mz, maxDiff = max_diff)
}
}

#' Fixed-width binning grouper
#'
#' Returns a grouper closure that assigns each m/z to a bin of width
#' `bin_width` (Da). Two values share a group iff `floor(mz / bin_width)`
#' is equal — fast, deterministic, and oblivious to clustering.
#'
#' Be aware of the bin-boundary artifact: values straddling a boundary
#' get split even if they're closer to each other than `bin_width`. For
#' tolerant grouping use [ppm_upgma()] / [abs_upgma()] instead.
#'
#' @param bin_width Numeric. Bin width in Da (typical: 0.01, 0.05, 0.1).
#' @return A `function(mz, target_mz)` satisfying the grouper contract.
#' @export
#' @examples
#' g <- bin_grouper(bin_width = 0.01)
#' g(c(100.001, 100.009, 100.015), target_mz = NA)
bin_grouper <- function(bin_width = 0.01) {
force(bin_width)
function(mz, target_mz = NULL) {
    as.integer(floor(mz / bin_width)) + 1L
}
}


#' Sanity-check a custom grouper against the contract
#'
#' Calls `grouper(mz_example, target_mz)` and asserts that the result is the
#' right length and has no more unique values than the input. Use this once
#' when wiring up a custom grouper; do not call inside hot loops.
#'
#' @param grouper A function with signature `function(mz, target_mz)`.
#' @param mz_example Numeric vector to probe with.
#' @param target_mz Numeric scalar passed to the grouper.
#' @return Invisibly `TRUE` when valid; throws otherwise.
#' @export
validate_grouper <- function(grouper,
                            mz_example = c(100.000, 100.001, 100.500, 200.000),
                            target_mz  = 100) {
if (!is.function(grouper)) stop("`grouper` must be a function.")
out <- grouper(mz_example, target_mz = target_mz)
if (length(out) != length(mz_example)) {
    stop(sprintf("grouper returned length %d, expected %d.",
                length(out), length(mz_example)))
}
if (length(unique(out)) > length(unique(mz_example))) {
    stop("grouper produced more unique groups than input m/z values.")
}
invisible(TRUE)
}