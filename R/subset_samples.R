#' Pick representative samples with apricot MaxCoverage
#'
#' A sample "covers" a feature when its intensity is at or above that
#' feature's `q_prob` quantile across samples (missing values count as
#' 0). apricot's `MaxCoverageSelection` then picks the `subset_size`
#' samples that together cover the most features. Needs Python with
#' `apricot-select` and `numpy` (via reticulate).
#'
#' @param msobject An `XcmsExperiment`/`XCMSnExp`.
#' @param msassay Samples × features intensity matrix;
#'   `rownames(msassay)` must be `basename(xcms::fileNames(msobject))`.
#' @param subset_size Number of samples to keep.
#' @param q_prob Per-feature quantile that defines "covered".
#' @param subsample If `FALSE`, return every sample and skip Python.
#' @return List with `msobject` (subsetted), `assay` (selected rows),
#'   `peaks` and `all_files` (both from the **unsubsetted** object, so
#'   `peakidx` and sample indices stay valid).
#' @keywords internal
#' @export
subset_samples_xcms <- function(
    msobject,
    msassay,
    subset_size = 100L,
    q_prob = 0.90,
    subsample = TRUE
) {
    if (!subsample) {
        return(list(
            msobject = msobject,
            assay    = msassay,
            peaks    = xcms::chromPeaks(msobject),
            all_files = xcms::fileNames(msobject)
        ))
    }

    if (!reticulate::py_module_available("apricot"))
        cli::cli_abort(c(
            "Python module {.pkg apricot} is not available to reticulate.",
            "i" = "reticulate is using {.file {reticulate::py_config()$python}}.",
            "i" = "The container sets RETICULATE_PYTHON to its own Python with apricot-select + numpy; a host run needs RETICULATE_PYTHON pointing at one."
        ))

    msassay[is.na(msassay)] <- 0
    thresholds <- apply(msassay, 2, quantile, probs = q_prob)
    A <- sweep(msassay, 2, thresholds, ">=")  # logical matrix, same dims as X
    A_np <- .pkg_py$np$array(A, dtype = "int32")

    selector <- .pkg_py$apricot$MaxCoverageSelection(n_samples = subset_size)  # note integer, not double
    selector$fit(A_np)
    selected_idx <- selector$ranking + 1L  # convert to 1-indexed for R
    assay_sub <- msassay[selected_idx,]
    peaks <- xcms::chromPeaks(msobject)

    object <- msobject[basename(xcms::fileNames(msobject)) %in% rownames(assay_sub)]
    return(list("msobject" = object, "assay" = assay_sub, "peaks" = peaks, all_files = xcms::fileNames(msobject)))
}

