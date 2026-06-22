#' Subset an XCMSnExp to a representative set of samples (apricot MaxCoverage)
#'
#' Picks `subset_size` files from `msobject` that collectively cover the
#' largest number of features present above the per-feature `q_prob`
#' quantile, using apricot's MaxCoverageSelection via reticulate. The full
#' `chromPeaks` table is returned alongside the subsetted msobject because
#' downstream `getBackendPeaks()` needs the original peak-row indexing
#' (`peakidx` from the parent featureDefinitions).
#'
#' @param msobject An `XCMSnExp`.
#' @param msassay  Numeric matrix, samples (rows) x features (cols).
#'   `rownames(msassay)` must match `basename(xcms::fileNames(msobject))`.
#' @param subset_size Integer. Number of samples to keep.
#' @param q_prob Numeric in (0, 1). Per-feature presence-threshold quantile.
#' @return A list with elements `msobject`, `assay`, `peaks`.
#' @export
subset_samples_xcms <- function(
    msobject,
    msassay,
    subset_size = 100L,
    q_prob = 0.90
) {
    msassay[is.na(msassay)] <- 0
    thresholds <- apply(msassay, 2, quantile, probs = q_prob)
    A <- sweep(msassay, 2, thresholds, ">=")  # logical matrix, same dims as X
    A_np <- .pkg_py$np$array(A, dtype = "int32")

    selector <- .pkg_py$apricot$MaxCoverageSelection(n_samples = subset_size)  # note integer, not double
    selector$fit(A_np)
    selected_idx <- selector$ranking + 1L  # convert to 1-indexed for R
    assay_sub <- msassay[selected_idx,]
    peaks <- xcms::chromPeaks(msobject) #|> as_tibble()

    object <- msobject[basename(xcms::fileNames(msobject)) %in% rownames(assay_sub)] 
    return(list("msobject" = object, "assay" = assay_sub, "peaks" = peaks))
}

