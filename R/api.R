

#' Process a single mzML file into a pseudo-spectra MGF
#'
#' @param path mzML file path.
#' @param tmp_dir working directory for intermediate parquet + outputs.
#' @param rt_range retention-time window (seconds) to retain.
#' @return cleaned pseudo-spectra data.table (invisibly written to MGF too).
#' @export
process_mzml <- function(path, tmp_dir, rt_range = c(30, 1200)){
    
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
    cleaned_spectra <- clean_pseudo_spectra(
        targets        = rownames(fd),
        file_directory = tmp_dir,
        fd             = as.data.frame(fd),
        ms1_ppm        = 10,
        grouper        = ppm_upgma(ppm = 30),
        buffer         = 50L,
        fdr_thresh     = 0.1,
        cor_thresh     = 0.1,
        min_count      = 10L,
        report         = FALSE,
        mgf_dir        = tmp_dir,
        progress       = TRUE
        )

    # --- export ---
    cli::cli_progress_step("Writing MGF")
    mgf_path <- file.path(tmp_dir, "all_features.mgf")
    export_mgf(cleaned_spectra, fd = fd, mgf_path = mgf_path, ion_mode = "POSITIVE")

    cli::cli_progress_done()
    cli::cli_alert_success("Finished {.file {basename(path)}}: {nrow(cleaned_spectra)} fragments across {length(unique(cleaned_spectra$feature_id))} features -> {.file {mgf_path}}")
    return(cleaned_spectra)
}


# ---- process_xcms ------------------------------------------------------------

#' Build pseudo-spectra MGF from a grouped XcmsExperiment
#'
#' One user-facing entry point for the xcms → MGF pipeline. Runs
#' end-to-end for interactive use (`stage = "all"`), or as a single
#' stage (`prepare` / `sample` / `finalize`) so that a DAG can
#' parallelize per-sample work.
#'
#' # On-disk contract
#'
#' All stages share a `tmp_dir` with a fixed layout — the DAG can
#' declare dependencies on these paths without asking R what was written:
#'
#' \itemize{
#'   \item `<tmp_dir>/manifest.csv` — `sample_idx,filename` (written by
#'     `prepare`; used by `sample` to resolve `sample_file`).
#'   \item `<tmp_dir>/mssubset.rds` — the [subset_samples_xcms()] output.
#'   \item `<tmp_dir>/fd.rds` — the `featureDefinitions()` frame.
#'   \item `<tmp_dir>/sample_<idx>_pool.parquet`,
#'         `<tmp_dir>/sample_<idx>_pointers.parquet` — one pair per
#'     sample (written by [Ms2PseudoRaw()] via [build_parquets_xcms()]).
#' }
#'
#' Split-stage runs must share the same `tmp_dir`.
#'
#' @param input An `XcmsExperiment`/`XCMSnExp`, or a path to one saved
#'   as `.rds` / `.qs` / `.rda`. Must have gone through
#'   `xcms::groupChromPeaks()`.
#' @param output MGF path. Required for `stage = "all"` or `"finalize"`.
#' @param tmp_dir Working directory for the on-disk contract. `NULL`
#'   only allowed for `stage = "all"` (creates + cleans a `tempfile`).
#'   Required otherwise.
#' @param stage One of `"all"` (default, interactive), `"prepare"`,
#'   `"sample"`, `"finalize"`.
#' @param sample_idx Integer index into `manifest$sample_idx` for
#'   `stage = "sample"`. Takes precedence over `sample_file`.
#' @param sample_file Basename of an mzML for `stage = "sample"`;
#'   resolved to an idx via the manifest.
#' @param subsample Logical. When `TRUE`, apricot picks `subset_size`
#'   representative samples; `FALSE` (default) uses every sample.
#' @param subset_size,q_prob Passed to [subset_samples_xcms()] when
#'   `subsample = TRUE`.
#' @param top_n Integer. Top-N samples per feature to build parquets
#'   for (passed to [build_parquets_xcms()]).
#' @param ppm Fragment/precursor tolerance in ppm.
#' @param rt_buffer Seconds added to each side of every peak's RT
#'   window (peak-alignment slack).
#' @param grouper Fragment-grouping function; `NULL` uses
#'   `ppm_upgma(ppm = ppm)`.
#' @param ion_mode `"POSITIVE"` or `"NEGATIVE"`; written to MGF headers.
#' @param ... Forwarded to [clean_pseudo_spectra()] (e.g. `fdr_thresh`,
#'   `cor_thresh`, `min_count`, `buffer`).
#' @return Invisibly the `tmp_dir` path (useful for chaining stages).
#' @family api
#' @export
process_xcms <- function(
    input,
    output     = NULL,
    tmp_dir    = NULL,
    stage      = c("all", "prepare", "sample", "finalize"),
    sample_idx = NULL,
    sample_file = NULL,
    subsample  = FALSE,
    subset_size = 100L,
    q_prob     = 0.90,
    top_n      = 10L,
    ppm        = 20,
    rt_buffer  = 5,
    grouper    = NULL,
    ion_mode   = "POSITIVE",
    ...
) {
    stage <- match.arg(stage)
    if (is.null(grouper)) grouper <- ppm_upgma(ppm = ppm)

    if (is.null(tmp_dir)) {
        if (stage != "all")
            cli::cli_abort("{.arg tmp_dir} is required when {.code stage != \"all\"}.")
        tmp_dir <- tempfile("carpediams_")
        dir.create(tmp_dir, recursive = TRUE)
        on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    } else {
        dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    }

    if (stage %in% c("all", "prepare"))
        .prepare_stage(input, tmp_dir, subsample, subset_size, q_prob)

    if (stage == "sample") {
        sid <- .resolve_sample_id(sample_idx, sample_file, tmp_dir)
        .sample_stage(tmp_dir, sid, top_n, ppm, rt_buffer)
    }

    if (stage == "all") {
        manifest <- utils::read.csv(file.path(tmp_dir, "manifest.csv"))
        cli::cli_progress_bar("Building per-sample parquets",
                              total = nrow(manifest))
        for (sid in manifest$sample_idx) {
            suppressMessages(.sample_stage(tmp_dir, sid, top_n, ppm, rt_buffer))
            cli::cli_progress_update()
        }
        cli::cli_progress_done()
    }

    if (stage %in% c("all", "finalize")) {
        if (is.null(output))
            cli::cli_abort("{.arg output} is required for stage {.val {stage}}.")
        .finalize_stage(tmp_dir, output, ppm, grouper, ion_mode, ...)
    }

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
.prepare_stage <- function(input, tmp_dir, subsample, subset_size, q_prob) {
    cli::cli_progress_step("Loading xcms object")
    xobj <- .load_xcms_input(input)
    fd <- xcms::featureDefinitions(xobj)
    if (nrow(fd) == 0)
        cli::cli_abort("Object has no features - run {.fn xcms::groupChromPeaks} first.")

    cli::cli_progress_step("Building assay matrix ({nrow(fd)} features)")
    msassay <- t(xcms::featureValues(xobj, method = "medret", value = "into"))
    rownames(msassay) <- basename(xcms::fileNames(xobj))

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

    manifest <- data.frame(
        sample_idx = seq_along(xcms::fileNames(mssubset$msobject)),
        filename   = xcms::fileNames(mssubset$msobject),
        stringsAsFactors = FALSE
    )
    utils::write.csv(manifest, file.path(tmp_dir, "manifest.csv"),
                     row.names = FALSE)
    saveRDS(mssubset, file.path(tmp_dir, "mssubset.rds"))
    saveRDS(as.data.frame(fd), file.path(tmp_dir, "fd.rds"))

    cli::cli_alert_success(
        "Prepare complete: {nrow(manifest)} sample{?s}, {nrow(fd)} feature{?s} -> {.path {tmp_dir}}")
}

#' @keywords internal
#' @noRd
.sample_stage <- function(tmp_dir, sample_idx, top_n, ppm, rt_buffer) {
    mssubset <- readRDS(file.path(tmp_dir, "mssubset.rds"))
    fd       <- readRDS(file.path(tmp_dir, "fd.rds"))
    build_parquets_xcms(
        mssubset  = mssubset, fd = fd, out_dir = tmp_dir,
        sample_id = sample_idx,
        top_n     = top_n,
        ppm       = ppm,
        rt_buffer = rt_buffer
    )
}

#' @keywords internal
#' @noRd
.finalize_stage <- function(tmp_dir, output, ppm, grouper, ion_mode, ...) {
    fd <- readRDS(file.path(tmp_dir, "fd.rds"))
    cli::cli_progress_step("Filtering pseudo-spectra ({nrow(fd)} features)")
    cleaned <- clean_pseudo_spectra(
        targets        = rownames(fd),
        file_directory = tmp_dir,
        fd             = fd,
        ms1_ppm        = ppm,
        grouper        = grouper,
        report         = FALSE,
        mgf_dir        = NULL,
        progress       = FALSE,
        ...
    )

    cli::cli_progress_step("Writing MGF to {.file {output}}")
    dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
    export_mgf(cleaned, fd = fd, mgf_path = output, ion_mode = ion_mode)

    cli::cli_progress_done()
    cli::cli_alert_success(
        "Done: {nrow(cleaned)} fragments across {length(unique(cleaned$feature_id))} feature{?s} -> {.file {output}}")
    invisible(cleaned)
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
