#' Consensus m/z contract
#'
#' All consensus-m/z closures in this package share one signature so they
#' are drop-in replaceable in [process_pseudo_spec()]:
#'
#' \preformatted{
#'   function(mz, int) -> single numeric consensus m/z
#' }
#'
#' Called once per fragment group, after grouping; the result becomes the
#' group's `fragment_mz`. Should return `NA_real_`, not error, when the
#' group has no finite values.
#'
#' @name consensus-mz-contract
NULL

#' Intensity-weighted FWHM-core consensus m/z (the package default)
#'
#' Intensity-weighted mean m/z of the centroids with intensity
#' `>= frac * max(intensity)` in the group, with weights capped at
#' `Isat`. Low-intensity centroids are dropped because their m/z is
#' least precise. Falls back to the unweighted mean of those centroids
#' if all weights are zero.
#'
#' @param frac Fraction of the group's maximum intensity a centroid
#'   needs to be included (`0.5` = half maximum).
#' @param Isat Weight cap, so a few saturated centroids cannot dominate.
#' @return A consensus-mz closure - see [consensus-mz-contract].
#' @seealso [process_pseudo_spec()]
#' @export
#' @examples
#' f <- consensus_mz_weighted(frac = 0.5, Isat = 1e4)
#' f(c(100.000, 100.001, 100.500), c(1000, 800, 50))
consensus_mz_weighted <- function(frac = 0.5, Isat = 1e4) {
  force(frac); force(Isat)
  function(mz, int) .consensus_mz(mz, int, Isat = Isat, frac = frac)
}
