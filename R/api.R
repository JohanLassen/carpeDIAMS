

#' Build a pseudo-spectra MGF from one mzML file
#'
#' Single-sample route with no prior feature list. Steps:
#'
#' 1. Read the file and keep spectra inside `rt_range`.
#' 2. If no spectrum is labelled MS2, assign MS levels by scan-index
#'    parity (bbCID alternates MS1/MS2). Aborts if the split is badly
#'    unbalanced.
#' 3. Detect chromatographic peaks with centWave (fixed parameters:
#'    `peakwidth = c(2, 30)`, `snthresh = 6`, `ppm = 12`,
#'    `prefilter = c(3, 1000)`, `mzdiff = 0.01`).
#' 4. Write pool/pointer parquets to `tmp_dir` and build one
#'    pseudo-spectrum per peak with [clean_pseudo_spectra()].
#' 5. Write the MGF with [export_mgf()].
#'
#' @param path mzML file path.
#' @param tmp_dir Directory for the intermediate parquets.
#' @param output MGF path. Required: the function stops before reading
#'   any data if it is `NULL`.
#' @param rt_range Retention-time window (seconds) to keep.
#' @param grouper A grouper closure (see [grouper-contract]); `NULL`
#'   uses `ppm_upgma(ppm = 30)`.
#' @param group_statistics_fn A statistics closure (see [stats-contract]);
#'   `NULL` uses `correlation_statistics(min_count = 10)`.
#' @param filter_fn A filter closure (see [filter-contract]); `NULL`
#'   uses `correlation_filtering(cor_thresh = 0.1, fdr_thresh = 0.1)`.
#' @return The [clean_pseudo_spectra()] table. The MGF is written as a
#'   side effect.
#' @family api
#' @export
process_mzml <- function(path,
                         tmp_dir,
                         output              = NULL,
                         rt_range            = c(30, 1200),
                         grouper             = NULL,
                         group_statistics_fn = NULL,
                         filter_fn           = NULL){
    
    if (is.null(output)) cli::cli_abort("No output provided.")

    cli::cli_h1("Processing {.file {basename(path)}}")
    cli::cli_progress_step("Loading spectra and filtering RT to {rt_range[1]}-{rt_range[2]}s")

    sps <- Spectra(path, backend = MsBackendMzR()) |> Spectra::filterRt(rt_range)
    if (length(sps) == 0) cli::cli_abort("No spectra in RT window.")

    # --- bbCID: assign MS2 by frame structure if not already annotated ---
    if (!2 %in% Spectra::msLevel(sps)) {
        cli::cli_alert_warning("No MS2 annotation found - assigning by frame parity (bbCID).")
        new_lvl <- as.integer((Spectra::scanIndex(sps) %% 2) + 1L)
        # sanity: expect a roughly balanced split, else the parity assumption is wrong
        tab <- table(new_lvl)
        if (length(tab) < 2 || min(tab) / max(tab) < 0.3)
        cli::cli_abort("Frame-parity MS-level assignment looks wrong (split: {paste(tab, collapse='/')}). Check acquisition / collision energy.")
        sps$msLevel <- new_lvl   # accessor, not slot-poking
        cli::cli_alert_info("Assigned MS1/MS2 by parity: {tab[1]} / {tab[2]} scans.")
    }
    
    cli::cli_progress_step("Detecting chromatographic peaks (centWave)")
    cwp <- CentWaveParam(
        peakwidth = c(2,30),
        snthresh  = 6,
        ppm       = 12,
        prefilter = c(3, 1000),
        mzdiff    = 0.01)

    mse <- MsExperiment()
    spectra(mse) <- sps

    sampleData(mse) <- DataFrame(sample_name = unique(dataOrigin(sps)))   # MUST equal dataOrigin
    mse <- linkSampleData(mse, with = "sampleData.sample_name = spectra.dataOrigin")
    sps <- findChromPeaks(mse, cwp)
    rm(mse);gc()

    fd <- chromPeaks(sps)
    rownames(fd) <- rownames(chromPeakData(sps))
    if (nrow(fd) == 0) cli::cli_abort("centWave found no peaks - check parameters.")
    cli::cli_alert_info("Found {nrow(fd)} chromatographic peaks.")
    
    spectra_df <- as.data.frame(spectraData(spectra(sps))) 
    fd_split <- fd |> as.data.frame() |> split(seq_len(nrow(fd)))
    names(fd_split) <- rownames(fd)

    cli::cli_progress_step("Extracting scan info and building pseudo-raw pool")
    scan_info <- extractScanInfo(feature_peaks = fd_split, filenames = fileNames(sps), spectra_df = spectra_df)
    spectra <- Ms2PseudoRaw(
        single_sample = scan_info,
        spectra_object = spectra(sps),
        maxDiff = 0.005,
        output_dir = tmp_dir,
        return_spectra = TRUE
    )
    
    cli::cli_progress_step("Building and filtering pseudo-spectra ({nrow(fd)} features)")
    if (is.null(grouper))             grouper             <- ppm_upgma(ppm = 30)
    if (is.null(group_statistics_fn)) group_statistics_fn <- correlation_statistics(min_count = 10L)
    if (is.null(filter_fn))           filter_fn           <- correlation_filtering(cor_thresh = 0.1, fdr_thresh = 0.1)

    cleaned_spectra <- clean_pseudo_spectra(
        targets             = rownames(fd),
        file_directory      = tmp_dir,
        fd                  = as.data.frame(fd),
        grouper             = grouper,
        group_statistics_fn = group_statistics_fn,
        filter_fn           = filter_fn,
        buffer              = 50L,
        report              = FALSE,
        mgf_dir             = NULL,
        progress            = TRUE
        )

    # --- export ---
    mgf_path <- if (is.null(output)) file.path(tmp_dir, "all_features.mgf") else output
    dir.create(dirname(mgf_path), showWarnings = FALSE, recursive = TRUE)
    cli::cli_progress_step("Writing MGF to {.file {mgf_path}}")
    export_mgf(cleaned_spectra, fd = fd, mgf_path = mgf_path, ion_mode = "POSITIVE")

    cli::cli_progress_done()
    cli::cli_alert_success("Finished {.file {basename(path)}}: {nrow(cleaned_spectra)} fragments across {length(unique(cleaned_spectra$feature_id))} features -> {.file {mgf_path}}")
    return(cleaned_spectra)
}


# ---- process_xcms ------------------------------------------------------------

#' Build a pseudo-spectra MGF from a grouped xcms object
#'
#' Multi-sample route. Every feature in `featureDefinitions()` gets one
#' MGF block, built from its top-`top_n` most intense samples. Run it in
#' one call (`stage = "all"`) or one stage at a time so a scheduler can
#' parallelize the per-sample and per-batch work.
#'
#' # Stages
#'
#' \describe{
#'   \item{`prepare`}{Loads `input`, builds the feature × sample
#'     intensity matrix (`featureValues(method = "medret", value =
#'     "into")`), optionally subsamples with apricot, and ranks the
#'     top-`top_n` samples per feature (ties can add more). Writes
#'     `manifest.csv`, `mssubset.rds`, `fd.rds`, `sample_topfeatures.rds`.}
#'   \item{`sample`}{For one sample: extracts the MS1/MS2 scans inside
#'     each qualifying feature's chromatographic peak (± `rt_buffer`)
#'     and writes `sample_<idx>_pool.parquet` (MS2) and
#'     `sample_<idx>_pointers.parquet` (MS1 trace per feature).}
#'   \item{`filter`}{For one contiguous slice (`batch_idx` of
#'     `n_batches`) of the features: groups, scores and filters
#'     fragments with [clean_pseudo_spectra()] and writes
#'     `features_batch_<batch_idx>.parquet`.}
#'   \item{`finalize`}{If any `features_batch_*.parquet` exist in
#'     `tmp_dir`, binds **all** of them; otherwise runs the filtering for
#'     every feature in-process. Then writes the MGF with [export_mgf()].}
#'   \item{`all`}{`prepare`, every `sample`, then `finalize` (in-process
#'     filtering, no batches).}
#' }
#'
#' All stages of one run must share `tmp_dir`.
#'
#' @param input An `XcmsExperiment`/`XCMSnExp` (after
#'   `xcms::groupChromPeaks()`), or a `.rds` / `.qs` / `.rda` path to
#'   one. Only read by `prepare`.
#' @param output MGF path. Required for `"all"` and `"finalize"`.
#' @param tmp_dir Working directory. `NULL` is only allowed for
#'   `stage = "all"`, which then uses a `tempfile()` deleted on exit.
#' @param stage One of `"all"`, `"prepare"`, `"sample"`, `"filter"`,
#'   `"finalize"`.
#' @param sample_idx Row of `manifest.csv` to process at `"sample"`.
#'   Takes precedence over `sample_file`.
#' @param sample_file mzML basename to process at `"sample"`; looked up
#'   in `manifest.csv`.
#' @param subsample If `TRUE`, keep `subset_size` samples chosen by
#'   apricot `MaxCoverageSelection` (needs Python with `apricot-select`
#'   and `numpy`). See [subset_samples_xcms()].
#' @param subset_size,q_prob Passed to [subset_samples_xcms()].
#' @param top_n Number of most intense samples per feature to extract.
#'   Used by `prepare` only.
#' @param ppm Used twice: at `sample`, the MS1 m/z tolerance for the
#'   precursor trace (converted to Da at the sample's largest feature
#'   m/z); and for the default grouper `ppm_upgma(ppm = ppm)`.
#' @param rt_buffer Seconds added to both sides of each peak's
#'   `rtmin`/`rtmax` at `sample`.
#' @param grouper A grouper closure (see [grouper-contract]); `NULL`
#'   uses `ppm_upgma(ppm = ppm)`.
#' @param group_statistics_fn A statistics closure (see [stats-contract]);
#'   `NULL` uses `correlation_statistics(min_count = min_count)`.
#' @param filter_fn A filter closure (see [filter-contract]); `NULL`
#'   uses `correlation_filtering(cor_thresh, fdr_thresh)`.
#' @param fdr_thresh,cor_thresh,min_count Only used to build the default
#'   `filter_fn` / `group_statistics_fn`.
#' @param batch_idx,n_batches For `"filter"`: process slice `batch_idx`
#'   (1-based) of `n_batches`. A slice past the last feature writes an
#'   empty parquet.
#' @param ion_mode `"POSITIVE"` or `"NEGATIVE"`, written as `IONMODE=`.
#' @param ... Forwarded to [clean_pseudo_spectra()] at `"filter"` and
#'   `"finalize"` (e.g. `buffer`, `consensus_mz_fn`, `apex_source`,
#'   `window_mode`).
#' @return Invisibly `tmp_dir`.
#' @family api
#' @export
process_xcms <- function(
    input,
    output              = NULL,
    tmp_dir             = NULL,
    stage               = c("all", "prepare", "sample", "filter", "finalize"),
    sample_idx          = NULL,
    sample_file         = NULL,
    subsample           = FALSE,
    subset_size         = 100L,
    q_prob              = 0.90,
    top_n               = 10L,
    ppm                 = 20,
    rt_buffer           = 5,
    grouper             = NULL,
    group_statistics_fn = NULL,
    filter_fn           = NULL,
    fdr_thresh          = 0.1,
    cor_thresh          = 0.1,
    min_count           = 10L,
    batch_idx           = NULL,
    n_batches           = 1L,
    ion_mode            = "POSITIVE",
    ...
) {
    stage <- match.arg(stage)
    .check_stage_args(stage, output, tmp_dir, sample_idx, sample_file, batch_idx)
    if (is.null(grouper))             grouper             <- ppm_upgma(ppm = ppm)
    if (is.null(group_statistics_fn)) group_statistics_fn <- correlation_statistics(min_count = min_count)
    if (is.null(filter_fn))           filter_fn           <- correlation_filtering(cor_thresh = cor_thresh,
                                                                                    fdr_thresh = fdr_thresh)

    if (is.null(tmp_dir)) {
        .warn_temporary_tmp_dir()
        tmp_dir <- tempfile("carpediams_")
        dir.create(tmp_dir, recursive = TRUE)
        on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    } else {
        dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    }

    if (stage %in% c("all", "prepare"))
        .prepare_stage(input, tmp_dir, subsample, subset_size, q_prob, top_n)

    if (stage == "sample") {
        sid <- .resolve_sample_id(sample_idx, sample_file, tmp_dir)
        .sample_stage(tmp_dir, sid, ppm, rt_buffer)
    }

    if (stage == "all") {
        manifest <- utils::read.csv(file.path(tmp_dir, "manifest.csv"))
        cli::cli_progress_bar("Building per-sample parquets",
                              total = nrow(manifest))
        for (sid in manifest$sample_idx) {
            suppressMessages(.sample_stage(tmp_dir, sid, ppm, rt_buffer))
            cli::cli_progress_update()
        }
        cli::cli_progress_done()
    }

    if (stage == "filter")
        .filter_stage(tmp_dir, batch_idx, n_batches, grouper,
                      group_statistics_fn, filter_fn, ...)

    if (stage %in% c("all", "finalize"))
        .finalize_stage(tmp_dir, output, ppm, grouper,
                        group_statistics_fn, filter_fn,
                        ion_mode, ...)

    invisible(tmp_dir)
}


# ---- internal stage helpers --------------------------------------------------

#' @keywords internal
#' @noRd
.load_xcms_input <- function(input) {
    if (inherits(input, c("XcmsExperiment", "XCMSnExp"))) return(input)
    if (!is.character(input) || length(input) != 1 || !file.exists(input))
        cli::cli_abort("Input must be an {.cls XcmsExperiment} or a path (.rds/.qs/.rda).")
    ext <- tolower(tools::file_ext(input))
    obj <- switch(ext,
        rds = readRDS(input),
        qs  = {
            if (!requireNamespace("qs", quietly = TRUE))
                cli::cli_abort("Install {.pkg qs} to read .qs files.")
            qs::qread(input)
        },
        rda = {
            e <- new.env(); load(input, envir = e)
            get(ls(e)[1], envir = e)
        },
        cli::cli_abort("Unsupported extension {.val {ext}} - use .rds/.qs/.rda.")
    )
    if (!inherits(obj, c("XcmsExperiment", "XCMSnExp")))
        cli::cli_abort("Loaded object is {.cls {class(obj)[1]}}, not an XcmsExperiment.")
    obj
}

#' @keywords internal
#' @noRd
.prepare_stage <- function(input, tmp_dir, subsample, subset_size, q_prob, top_n) {
    cli::cli_progress_step("Loading xcms object")
    xobj  <- .load_xcms_input(input)
    fd_df <- as.data.frame(xcms::featureDefinitions(xobj))
    if (nrow(fd_df) == 0)
        cli::cli_abort("Object has no features - run {.fn xcms::groupChromPeaks} first.")

    cli::cli_progress_step("Building assay matrix ({nrow(fd_df)} features)")
    msassay <- t(xcms::featureValues(xobj, method = "medret", value = "into"))
    rownames(msassay) <- basename(xcms::fileNames(xobj))
    if (ncol(msassay) != nrow(fd_df))
        cli::cli_abort(
            "Assay matrix has {ncol(msassay)} feature column{?s} but
             featureDefinitions has {nrow(fd_df)} row{?s} - xcms object
             looks inconsistent.")

    # fd.rds keeps every feature, valid or not - a feature with no usable
    # precursor mz (e.g. an empty peakidx, or a hand-edited xcms object) is
    # not filtered out here. clean_pseudo_spectra() reports it per-target
    # as a DeconvERROR:MissingTargetInfo sentinel row at finalize(), so
    # every feature in fd unconditionally gets exactly one MGF block -
    # never silently dropped, and never carried as separate sidecar state.

    n_in <- nrow(msassay)
    cli::cli_progress_step(
        if (subsample) "Subsampling with apricot ({subset_size} of {n_in} samples)"
        else "Passing through {n_in} sample{?s} (no subsampling)"
    )
    mssubset <- subset_samples_xcms(
        xobj, msassay,
        subset_size = subset_size, q_prob = q_prob,
        subsample = subsample
    )

    # Top-N-per-feature ranking, computed once against mssubset$assay - i.e.
    # AFTER subset_samples_xcms() has resolved apricot-on (subsetted) or
    # apricot-off (full passthrough). Previously this exact computation was
    # redone from scratch inside every per-sample build_parquets_xcms() call;
    # relocating it here turns O(samples) full-matrix rescans into one.
    cli::cli_progress_step("Ranking top-{top_n} sample{?s} per feature")
    assay_tbl <- tibble::as_tibble(mssubset$assay, rownames = "name")
    colnames(assay_tbl) <- c("name", rownames(fd_df))
    topX <- assay_tbl |>
        tidyr::pivot_longer(cols = -"name", names_to = "feature_id") |>
        tidyr::replace_na(list(value = 0)) |>
        dplyr::group_by(feature_id) |>
        dplyr::slice_max(value, n = top_n) |>
        dplyr::ungroup()
    sample_topfeatures <- split(topX$feature_id, topX$name)
    rm(assay_tbl, topX)

    manifest <- data.frame(
        sample_idx = seq_along(xcms::fileNames(mssubset$msobject)),
        filename   = xcms::fileNames(mssubset$msobject),
        stringsAsFactors = FALSE
    )
    utils::write.csv(manifest, file.path(tmp_dir, "manifest.csv"),
                     row.names = FALSE)
    saveRDS(mssubset, file.path(tmp_dir, "mssubset.rds"))
    saveRDS(fd_df, file.path(tmp_dir, "fd.rds"))
    saveRDS(sample_topfeatures, file.path(tmp_dir, "sample_topfeatures.rds"))

    cli::cli_alert_success(
        "Prepare complete: {nrow(manifest)} sample{?s}, {nrow(fd_df)} feature{?s} -> {.path {tmp_dir}}")
}

#' @keywords internal
#' @noRd
.sample_stage <- function(tmp_dir, sample_idx, ppm, rt_buffer) {
    mssubset           <- readRDS(file.path(tmp_dir, "mssubset.rds"))
    fd                 <- readRDS(file.path(tmp_dir, "fd.rds"))
    sample_topfeatures <- readRDS(file.path(tmp_dir, "sample_topfeatures.rds"))

    sample_file         <- basename(xcms::fileNames(mssubset$msobject)[sample_idx])
    qualifying_features <- sample_topfeatures[[sample_file]]
    if (is.null(qualifying_features)) qualifying_features <- character(0)

    build_parquets_xcms(
        mssubset            = mssubset, fd = fd, out_dir = tmp_dir,
        sample_id           = sample_idx,
        qualifying_features = qualifying_features,
        ppm                 = ppm,
        rt_buffer           = rt_buffer
    )
}

#' @keywords internal
#' @noRd
.finalize_stage <- function(tmp_dir, output, ppm, grouper,
                            group_statistics_fn, filter_fn,
                            ion_mode, require_rt = FALSE, ...) {
    fd      <- readRDS(file.path(tmp_dir, "fd.rds"))
    batches <- list.files(tmp_dir, pattern = "^features_batch_.*\\.parquet$",
                          full.names = TRUE)

    dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)

    if (length(batches)) {
        cli::cli_progress_step("Combining {length(batches)} filter batch{?s}")
        cleaned <- data.table::rbindlist(
            lapply(batches, arrow::read_parquet),
            fill = TRUE
        )
    } else {
        cli::cli_progress_step("Filtering pseudo-spectra ({nrow(fd)} features)")
        cleaned <- clean_pseudo_spectra(
            targets             = rownames(fd),
            file_directory      = tmp_dir,
            fd                  = fd,
            grouper             = grouper,
            group_statistics_fn = group_statistics_fn,
            filter_fn           = filter_fn,
            require_rt          = require_rt,
            report              = FALSE,
            mgf_dir             = NULL,
            progress            = FALSE,
            ...
        )
    }

    cli::cli_progress_step("Writing MGF to {.file {output}}")
    export_mgf(cleaned, fd = fd, mgf_path = output, ion_mode = ion_mode)

    cli::cli_progress_done()
    cli::cli_alert_success(
        "Done: {nrow(fd)} feature{?s} -> {.file {output}} (see DeconvERROR titles for any without a real spectrum)")
    invisible(cleaned)
}

#' @keywords internal
#' @noRd
.filter_stage <- function(tmp_dir, batch_idx, n_batches, grouper,
                          group_statistics_fn, filter_fn, require_rt = FALSE, ...) {
    fd <- readRDS(file.path(tmp_dir, "fd.rds"))
    cli::cli_progress_step(
        "Filtering pseudo-spectra (batch {batch_idx}/{n_batches})")
    cleaned <- clean_pseudo_spectra(
        targets             = rownames(fd),
        file_directory      = tmp_dir,
        fd                  = fd,
        grouper             = grouper,
        group_statistics_fn = group_statistics_fn,
        filter_fn           = filter_fn,
        require_rt          = require_rt,
        batch_idx           = batch_idx,
        n_batches           = n_batches,
        report              = FALSE,
        mgf_dir             = NULL,
        progress            = TRUE,
        ...
    )
    out_path <- file.path(tmp_dir,
                          sprintf("features_batch_%d.parquet", as.integer(batch_idx)))
    arrow::write_parquet(cleaned, out_path)
    cli::cli_progress_done()
    cli::cli_alert_success(
        "Wrote batch {batch_idx}/{n_batches}: {nrow(cleaned)} row{?s} -> {.path {out_path}}")
    invisible(cleaned)
}

#' Abort before any work if `stage` is missing an argument it needs.
#' @keywords internal
#' @noRd
.check_stage_args <- function(stage, output, tmp_dir, sample_idx,
                              sample_file, batch_idx) {
    if (stage != "all" && is.null(tmp_dir))
        cli::cli_abort("{.arg tmp_dir} is required when {.code stage != \"all\"}. Nothing was run.")
    if (stage == "sample" && is.null(sample_idx) && is.null(sample_file))
        cli::cli_abort("Provide {.arg sample_idx} or {.arg sample_file} for {.code stage = \"sample\"}. Nothing was run.")
    if (stage == "filter" && is.null(batch_idx))
        cli::cli_abort("{.code stage = \"filter\"} requires {.arg batch_idx}. Nothing was run.")
    if (stage %in% c("all", "finalize") && is.null(output)) {
        has_results <- !is.null(tmp_dir) &&
            length(list.files(tmp_dir, pattern = "\\.parquet$")) > 0
        cli::cli_abort(c(
            "{.arg output} is required for stage {.val {stage}}. Nothing was run.",
            if (has_results) c("i" = "Results already in {.path {tmp_dir}} can be written to an MGF with {.code stage = \"finalize\"} and an {.arg output} path.")
        ))
    }
    invisible(TRUE)
}

#' @keywords internal
#' @noRd
.warn_temporary_tmp_dir <- function() {
    cli::cli_alert_warning(
        "No {.arg tmp_dir} given: intermediate results are deleted when this
         call ends, also if it fails. Pass {.arg tmp_dir} to keep them; a
         failed run can then be finished with {.code stage = \"finalize\"}.")
}

#' @keywords internal
#' @noRd
.resolve_sample_id <- function(sample_idx, sample_file, tmp_dir) {
    if (!is.null(sample_idx)) return(as.integer(sample_idx))
    if (is.null(sample_file))
        cli::cli_abort("Provide {.arg sample_idx} or {.arg sample_file} for {.code stage = \"sample\"}.")
    manifest <- utils::read.csv(file.path(tmp_dir, "manifest.csv"))
    hit <- manifest$sample_idx[basename(manifest$filename) == basename(sample_file)]
    if (length(hit) == 0)
        cli::cli_abort("{.arg sample_file} {.val {sample_file}} not found in manifest.")
    as.integer(hit[1])
}


# ---- process_targetlist ------------------------------------------------------

#' Build a pseudo-spectra MGF from mzML file(s) and a target list
#'
#' Targeted route: no peak detection. Every target is extracted from
#' every file, in a fixed RT window around its `rt`. Same stages and
#' `tmp_dir` layout as [process_xcms()], except that `prepare` writes
#' only `manifest.csv` and `fd.rds` (no subsampling, no top-N ranking).
#'
#' Target-list handling:
#'
#' * Needs `mz` and `rt` columns; `feature` is optional (ids
#'   `FT01`, `FT02`, … are generated when missing).
#' * If the largest `rt` is below 60 it is taken as minutes and
#'   converted to seconds.
#' * Rows with missing `mz` or `rt` are kept and reported as
#'   `DeconvERROR:MissingTargetInfo` in the MGF.
#'
#' At `sample`, files with no MS2 labels get MS levels by scan-index
#' parity (as in [process_mzml()]). Samples where no target falls within
#' the window write no parquets; `stage = "all"` reports them.
#'
#' @inheritParams process_xcms
#' @param input One mzML path, a directory of mzMLs (non-recursive), or
#'   a character vector of mzML paths. Only read by `prepare`.
#' @param targetlist CSV path or data.frame (see Details). Required for
#'   `"all"` and `"prepare"`.
#' @param ppm Used twice: at `sample`, the MS1 m/z tolerance for the
#'   precursor trace (converted to Da at the largest target m/z); and
#'   for the default grouper `ppm_upgma(ppm = ppm)`.
#' @param rt_width Total width (seconds) of the extraction window,
#'   centred on each target's `rt`.
#' @param ... Forwarded to [clean_pseudo_spectra()] at `"filter"` and
#'   `"finalize"`. `require_rt` is always `TRUE` on this route.
#' @return Invisibly `tmp_dir`.
#' @family api
#' @export
process_targetlist <- function(
    input,
    output              = NULL,
    targetlist          = NULL,
    tmp_dir             = NULL,
    stage               = c("all", "prepare", "sample", "filter", "finalize"),
    sample_idx          = NULL,
    sample_file         = NULL,
    ppm                 = 20,
    rt_width            = 15,
    grouper             = NULL,
    group_statistics_fn = NULL,
    filter_fn           = NULL,
    fdr_thresh          = 0.1,
    cor_thresh          = 0.1,
    min_count           = 10L,
    batch_idx           = NULL,
    n_batches           = 1L,
    ion_mode            = "POSITIVE",
    ...
) {
    stage <- match.arg(stage)
    .check_stage_args(stage, output, tmp_dir, sample_idx, sample_file, batch_idx)
    if (stage %in% c("all", "prepare") && is.null(targetlist))
        cli::cli_abort("{.arg targetlist} is required for stage {.val {stage}}. Nothing was run.")
    if (is.null(grouper))             grouper             <- ppm_upgma(ppm = ppm)
    if (is.null(group_statistics_fn)) group_statistics_fn <- correlation_statistics(min_count = min_count)
    if (is.null(filter_fn))           filter_fn           <- correlation_filtering(cor_thresh = cor_thresh,
                                                                                    fdr_thresh = fdr_thresh)

    if (is.null(tmp_dir)) {
        .warn_temporary_tmp_dir()
        tmp_dir <- tempfile("carpediams_")
        dir.create(tmp_dir, recursive = TRUE)
        on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    } else {
        dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    }

    if (stage %in% c("all", "prepare"))
        .prepare_targetlist_stage(input, targetlist, tmp_dir)

    if (stage == "sample") {
        sid <- .resolve_sample_id(sample_idx, sample_file, tmp_dir)
        .sample_targetlist_stage(tmp_dir, sid, ppm, rt_width)
    }

    if (stage == "all") {
        manifest <- utils::read.csv(file.path(tmp_dir, "manifest.csv"))
        cli::cli_progress_bar("Building per-sample parquets",
                              total = nrow(manifest))
        skipped <- integer(0)
        for (sid in manifest$sample_idx) {
            ok <- suppressMessages(
                .sample_targetlist_stage(tmp_dir, sid, ppm, rt_width))
            if (!isTRUE(ok)) skipped <- c(skipped, sid)
            cli::cli_progress_update()
        }
        cli::cli_progress_done()
        # suppressMessages() above keeps the progress bar readable, but it also
        # hides the per-sample "no targets in window" notice. Report it here.
        if (length(skipped))
            cli::cli_alert_warning(
                "{length(skipped)} of {nrow(manifest)} samples produced no
                 pointers (sample_idx {skipped}): no target fell within
                 {rt_width}s of an MS1 scan.")
    }

    if (stage == "filter")
        .filter_stage(tmp_dir, batch_idx, n_batches, grouper,
                      group_statistics_fn, filter_fn, require_rt = TRUE, ...)

    if (stage %in% c("all", "finalize"))
        .finalize_stage(tmp_dir, output, ppm, grouper,
                        group_statistics_fn, filter_fn,
                        ion_mode, require_rt = TRUE, ...)

    invisible(tmp_dir)
}


# ---- internal stage helpers (targetlist) -------------------------------------

#' @keywords internal
#' @noRd
.resolve_targetlist_input <- function(input) {
    if (length(input) > 1) {
        missing <- input[!file.exists(input)]
        if (length(missing))
            cli::cli_abort("Missing file{?s}: {.path {missing}}")
        return(input)
    }
    if (dir.exists(input)) {
        files <- list.files(input, pattern = "\\.mzML$",
                            full.names = TRUE, ignore.case = TRUE)
        if (!length(files))
            cli::cli_abort("No .mzML files in {.path {input}}")
        return(files)
    }
    if (file.exists(input)) return(input)
    cli::cli_abort("Input must be an mzML file, directory, or vector of mzML paths.")
}

#' @keywords internal
#' @noRd
.load_targetlist <- function(targetlist) {
    tl <- if (is.character(targetlist)) {
        readr::read_csv(targetlist, show_col_types = FALSE)
    } else as.data.frame(targetlist)
    miss <- setdiff(c("mz", "rt"), colnames(tl))
    if (length(miss)) cli::cli_abort("targetlist missing column{?s}: {miss}")

    # Assign ids to every row (including any with missing mz/rt) up front,
    # before the unit-heuristic below might otherwise be tempted to drop
    # anything - a row with missing mz/rt is kept and given a stable id
    # traceable to its original CSV row; clean_pseudo_spectra() reports it
    # as a DeconvERROR:MissingTargetInfo sentinel later, it is never
    # filtered out here.
    if (!"feature" %in% colnames(tl))
        tl$feature <- sprintf("FT%0*d", nchar(nrow(tl)), seq_len(nrow(tl)))
    tl$feature <- as.character(tl$feature)

    if (any(!is.na(tl$rt)) && max(tl$rt, na.rm = TRUE) < 60) {
        cli::cli_alert_info(
            "RT looks like minutes (max = {round(max(tl$rt, na.rm = TRUE), 2)}); converting to seconds.")
        tl$rt <- tl$rt * 60
    }
    tl
}

#' @keywords internal
#' @noRd
.prepare_targetlist_stage <- function(input, targetlist, tmp_dir) {
    cli::cli_progress_step("Resolving input files")
    files <- .resolve_targetlist_input(input)

    cli::cli_progress_step("Loading target list")
    tl <- .load_targetlist(targetlist)

    manifest <- data.frame(
        sample_idx = seq_along(files),
        filename   = files,
        stringsAsFactors = FALSE
    )
    utils::write.csv(manifest, file.path(tmp_dir, "manifest.csv"),
                     row.names = FALSE)

    fd <- data.frame(
        feature_id = tl$feature,
        mzmed      = tl$mz,
        rtmed      = tl$rt,
        row.names  = tl$feature,
        stringsAsFactors = FALSE
    )
    saveRDS(fd, file.path(tmp_dir, "fd.rds"))

    cli::cli_alert_success(
        "Prepare complete: {nrow(manifest)} file{?s}, {nrow(fd)} target{?s} -> {.path {tmp_dir}}")
}

#' @keywords internal
#' @noRd
.sample_targetlist_stage <- function(tmp_dir, sample_idx, ppm, rt_width) {
    manifest <- utils::read.csv(file.path(tmp_dir, "manifest.csv"))
    file_name <- manifest$filename[manifest$sample_idx == sample_idx]
    if (length(file_name) != 1)
        cli::cli_abort("sample_idx {.val {sample_idx}} not found in manifest.")

    fd <- readRDS(file.path(tmp_dir, "fd.rds"))
    if (nrow(fd) == 0) {
        # A targetlist can "prepare" successfully into zero usable targets
        # (e.g. every row's mz/rt was NA - see .load_targetlist()'s filter).
        # That is a clean, reportable outcome, not a crash: without this
        # guard the data.frame() call below fails with an opaque
        # "differing number of rows" error, because file_name/sample_idx
        # are length 1 but fd$feature_id etc. are length 0. Same reason
        # code clean_pseudo_spectra() uses for this condition via
        # .deconv_error_row() - see ?clean_pseudo_spectra.
        message("Sample ", sample_idx, ": MissingTargetInfo (target list has no usable mz/rt).")
        return(invisible("MissingTargetInfo"))
    }

    sps <- Spectra::Spectra(file_name, backend = Spectra::MsBackendMzR())
    if (!2 %in% Spectra::msLevel(sps)) {
        new_lvl <- as.integer((Spectra::scanIndex(sps) %% 2) + 1L)
        tab <- table(new_lvl)
        if (length(tab) < 2 || min(tab) / max(tab) < 0.3)
            cli::cli_abort(
                "Frame-parity MS-level assignment for {.file {basename(file_name)}} looks wrong (split: {paste(tab, collapse='/')}).")
        sps$msLevel <- new_lvl
    }

    spectra_df <- as.data.frame(Spectra::spectraData(sps)) |>
        dplyr::filter(msLevel == 1) |>
        dplyr::select(dataOrigin, rtime, acquisitionNum)

    scan_info <- data.frame(
        feature  = fd$feature_id,
        mz       = fd$mzmed,
        rt       = fd$rtmed,
        filename = file_name,
        sample   = sample_idx,
        stringsAsFactors = FALSE
    ) |>
        dplyr::inner_join(spectra_df,
                          by = c("filename" = "dataOrigin"),
                          relationship = "many-to-many") |>
        dplyr::filter(rtime >= rt - rt_width / 2,
                      rtime <= rt + rt_width / 2) |>
        dplyr::group_by(feature, sample, filename, mz, rt) |>
        dplyr::summarise(
            scan_indices    = list(acquisitionNum),
            apex_scan_index = acquisitionNum[which.min(abs(rtime - rt))],
            .groups = "drop"
        ) |>
        tidyr::unnest(scan_indices)

    if (nrow(scan_info) == 0) {
        message("Sample ", sample_idx, ": MissingMs1Signal (no targets fell in the ", rt_width, "s window).")
        return(invisible("MissingMs1Signal"))
    }

    max_diff <- max(fd$mzmed, na.rm = TRUE) * ppm * 1e-6
    res <- Ms2PseudoRaw(
        single_sample  = scan_info,
        spectra_object = sps,
        output_dir     = tmp_dir,
        maxDiff        = max_diff
    )
    # Ms2PseudoRaw() returns NULL without writing when the acquisition-number
    # window selects no spectra - same MissingMs1Signal condition, just
    # caught one level deeper (a target can pass the RT-window filter above
    # via one scan yet still have no peaksData in that scan).
    wrote <- file.exists(file.path(tmp_dir, paste0("sample_", sample_idx, "_pointers.parquet")))
    if (!wrote) {
        message("Sample ", sample_idx, ": MissingMs1Signal (Ms2PseudoRaw wrote no pointer rows).")
        return(invisible("MissingMs1Signal"))
    }
    invisible(TRUE)
}


# ---- carpediams (unified dispatch) -------------------------------------------

#' Build a pseudo-spectra MGF (dispatcher)
#'
#' Picks [process_xcms()], [process_targetlist()] or [process_mzml()]
#' from `input` and `targetlist`, and passes everything else on.
#'
#' # Routing
#'
#' | `input`                                         | `targetlist` | Route                   |
#' |:------------------------------------------------|:-------------|:------------------------|
#' | `XcmsExperiment`/`XCMSnExp`, or `.rds`/`.qs`/`.rda` path | ignored (warning) | [process_xcms()] |
#' | mzML file, mzML directory, or vector of paths   | given        | [process_targetlist()]  |
#' | one mzML file                                   | `NULL`       | [process_mzml()]        |
#' | mzML directory or vector of paths               | `NULL`       | error                   |
#'
#' # Pluggable steps
#'
#' * `grouper` ([grouper-contract]): [ppm_upgma()], [abs_upgma()],
#'   [iaac_grouper()], [pelt_grouper()], [bin_grouper()].
#' * `group_statistics_fn` ([stats-contract]): [correlation_statistics()],
#'   [apex_correlation_statistics()].
#' * `filter_fn` ([filter-contract]): [correlation_filtering()],
#'   [apex_correlation_filtering()] (needs [apex_correlation_statistics()]).
#' * `consensus_mz_fn` ([consensus-mz-contract]): [consensus_mz_weighted()].
#'   Pass it through `...`; xcms and target-list routes only.
#'
#' @param input Object or path; see Routing.
#' @param output MGF path. Not needed for the `"prepare"`, `"sample"`
#'   and `"filter"` stages.
#' @param targetlist `NULL`, or a CSV path / data.frame with `mz` and
#'   `rt` columns.
#' @param tmp_dir Working directory. `NULL` uses a `tempfile()` that is
#'   **deleted when the call returns**, so staged runs must pass a path.
#' @param grouper A grouper closure; `NULL` uses the route default
#'   (`ppm_upgma(ppm = 20)` for xcms/target-list via their `ppm`
#'   argument, `ppm_upgma(ppm = 30)` for mzML).
#' @param group_statistics_fn A statistics closure; `NULL` uses
#'   `correlation_statistics(min_count = 10)`.
#' @param filter_fn A filter closure; `NULL` uses
#'   `correlation_filtering(cor_thresh = 0.1, fdr_thresh = 0.1)`.
#' @param ... Passed to the routed function (`stage`, `ppm`, `top_n`,
#'   …). [process_mzml()] takes no extra arguments.
#' @return What the routed function returns: `tmp_dir` (invisibly) for
#'   the xcms and target-list routes, the pseudo-spectra table for mzML.
#' @family api
#' @examples
#' \dontrun{
#' # default — correlation-based filter
#' carpediams(input = "grouped.rds", output = "out.mgf")
#'
#' # swap the grouper
#' carpediams(input = "grouped.rds", output = "out.mgf",
#'            grouper = iaac_grouper(k = 1, p0 = 8, Isat = 1e5))
#'
#' # apex-agreement filter (needs the matching statistics)
#' carpediams(input = "grouped.rds", output = "out.mgf",
#'            group_statistics_fn = apex_correlation_statistics(),
#'            filter_fn           = apex_correlation_filtering(
#'                                    max_apex_shift_med = 0.5,
#'                                    max_apex_shift_iqr = 0.75))
#' }
#' @export
carpediams <- function(input, output,
                       targetlist          = NULL,
                       tmp_dir             = NULL,
                       grouper             = NULL,
                       group_statistics_fn = NULL,
                       filter_fn           = NULL,
                       ...) {
    if (missing(output)) output <- NULL
    route <- .detect_route(input, targetlist)
    cli::cli_h1("carpediams route: {.strong {route}}")

    if (route != "mzml" && is.null(tmp_dir)) {
        stage <- if ("stage" %in% ...names()) ...elt(match("stage", ...names())) else "all"
        if (!identical(stage, "all"))
            cli::cli_abort("{.arg tmp_dir} is required when {.code stage != \"all\"}. Nothing was run.")
        .warn_temporary_tmp_dir()
    }

    if (is.null(tmp_dir)) {
        tmp_dir <- tempfile("carpediams_")
        dir.create(tmp_dir, recursive = TRUE)
        on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    } else {
        dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    }

    switch(route,
        xcms       = process_xcms(
            input               = input,
            output              = output,
            tmp_dir             = tmp_dir,
            grouper             = grouper,
            group_statistics_fn = group_statistics_fn,
            filter_fn           = filter_fn,
            ...),
        targetlist = process_targetlist(
            input               = input,
            output              = output,
            targetlist          = targetlist,
            tmp_dir             = tmp_dir,
            grouper             = grouper,
            group_statistics_fn = group_statistics_fn,
            filter_fn           = filter_fn,
            ...),
        mzml       = process_mzml(
            path                = input,
            tmp_dir             = tmp_dir,
            output              = output,
            grouper             = grouper,
            group_statistics_fn = group_statistics_fn,
            filter_fn           = filter_fn,
            ...)
    )
}

#' @keywords internal
#' @noRd
.detect_route <- function(input, targetlist) {
    if (inherits(input, c("XcmsExperiment", "XCMSnExp"))) {
        if (!is.null(targetlist))
            cli::cli_alert_warning("{.arg targetlist} ignored for xcms object input.")
        return("xcms")
    }
    if (!is.character(input))
        cli::cli_abort("{.arg input} must be a character path/vector or an {.cls XcmsExperiment}.")

    if (length(input) > 1) {
        if (is.null(targetlist))
            cli::cli_abort("Multiple mzML paths require a {.arg targetlist} - or call {.fn process_mzml} per file.")
        return("targetlist")
    }
    if (dir.exists(input)) {
        if (is.null(targetlist))
            cli::cli_abort("A directory input requires a {.arg targetlist}.")
        return("targetlist")
    }
    if (!file.exists(input))
        cli::cli_abort("{.arg input} path does not exist: {.path {input}}")

    ext <- tolower(tools::file_ext(input))
    if (ext %in% c("rds", "qs", "rda")) {
        if (!is.null(targetlist))
            cli::cli_alert_warning("{.arg targetlist} ignored for serialized xcms input.")
        return("xcms")
    }
    if (ext == "mzml")
        return(if (is.null(targetlist)) "mzml" else "targetlist")

    cli::cli_abort("Unrecognized file extension {.val {ext}} for {.arg input}.")
}
