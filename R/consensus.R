#' Consensus m/z contract
#'
#' All consensus-m/z closures in this package share one signature so they
#' are drop-in replaceable in [process_pseudo_spec()]:
#'
#' \preformatted{
#'   function(mz, int) -> single numeric consensus m/z
#' }
#'
#' Called once per group at the pass-2 collapse step (per `(sample,
#' pass-1 group)`) and once per final group (writing `group_mz`). Must
#' return `NA_real_` for a group with no finite values, never error.
#'
#' @name consensus-mz-contract
#' @keywords internal
NULL

#' Intensity-weighted FWHM-core consensus m/z (the package default)
#'
#' Returns a consensus-mz closure wrapping the package's default
#' estimator: intensity-weighted mean over centroids at or above
#' `frac * max(intensity)` within the group, weight clipped at `Isat`.
#' See [process_pseudo_spec()] for why the FWHM core is used
#' (low-intensity centroids carry disproportionate m/z error).
#'
#' @param frac Numeric in (0, 1]. FWHM-core cutoff as a fraction of the
#'   group's max intensity.
#' @param Isat Numeric. Weight-clipping intensity (prevents a saturated
#'   apex from dominating the mean).
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
