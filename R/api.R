

#' Process a single mzML file into a pseudo-spectra MGF
#'
#' @param path mzML file path.
#' @param tmp_dir working directory for intermediate parquet + outputs.
#' @param output MGF path. `NULL` (default) writes to
#'   `file.path(tmp_dir, "all_features.mgf")` for backward compatibility.
#' @param rt_range retention-time window (seconds) to retain.
#' @return cleaned pseudo-spectra data.table (invisibly written to MGF too).
#' @family api
#' @export
process_mzml <- function(path, tmp_dir, output = NULL, rt_range = c(30, 1200)){
    
    cli::cli_h1("Processing {.file {basename(path)}}")
    cli::cli_progress_step("Loading spectra and filtering RT to {rt_range[1]}-{rt_range[2]}s")

    sps <- Spectra(path, backend = MsBackendMzR()) |> Spectra::filterRt(rt_range)
    if (length(sps) == 0) cli::cli_abort("No spectra in RT window.")
    if (output == NULL) cli::cli_abort("No output provided.")

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
        grouper        = ppm_upgma(ppm = 30),
        buffer         = 50L,
        fdr_thresh     = 0.1,
        cor_thresh     = 0.1,
        min_count      = 10L,
        report         = FALSE,
        mgf_dir        = NULL,
        progress       = TRUE
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
#'   \item `<tmp_dir>/features_batch_<i>.parquet` — one per `filter`
#'     batch (present only when the batched path is used; `finalize`
#'     auto-detects and aggregates them).
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
#'   `"sample"`, `"filter"`, `"finalize"`. `"filter"` is the batched
#'   scale-out path (see `batch_idx`); `"finalize"` auto-aggregates any
#'   batch parquets if present, else runs everything in-process.
#' @param sample_idx Integer index into `manifest$sample_idx` for
#'   `stage = "sample"`. Takes precedence over `sample_file`.
#' @param sample_file Basename of an mzML for `stage = "sample"`;
#'   resolved to an idx via the manifest.
#' @param batch_idx,n_batches For `stage = "filter"`: split `targets`
#'   into `n_batches` contiguous chunks and process chunk `batch_idx`
#'   (1-based). One `filter` DAG job per batch; `finalize` unions them.
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
#' @param fdr_thresh,cor_thresh,min_count Passed to
#'   [clean_pseudo_spectra()]. Defaults mirror [process_mzml()]
#'   (permissive) rather than [clean_pseudo_spectra()]'s stricter
#'   built-in defaults.
#' @param ion_mode `"POSITIVE"` or `"NEGATIVE"`; written to MGF headers.
#' @param ... Forwarded to [clean_pseudo_spectra()] (e.g. `buffer`).
#' @return Invisibly the `tmp_dir` path (useful for chaining stages).
#' @family api
#' @export
process_xcms <- function(
    input,
    output     = NULL,
    tmp_dir    = NULL,
    stage      = c("all", "prepare", "sample", "filter", "finalize"),
    sample_idx = NULL,
    sample_file = NULL,
    subsample  = FALSE,
    subset_size = 100L,
    q_prob     = 0.90,
    top_n      = 10L,
    ppm        = 20,
    rt_buffer  = 5,
    grouper    = NULL,
    fdr_thresh = 0.1,
    cor_thresh = 0.1,
    min_count  = 10L,
    batch_idx  = NULL,
    n_batches  = 1L,
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

    if (stage == "filter") {
        if (is.null(batch_idx))
            cli::cli_abort("stage = {.val filter} requires {.arg batch_idx}.")
        .filter_stage(tmp_dir, batch_idx, n_batches, grouper,
                      fdr_thresh, cor_thresh, min_count, ...)
    }

    if (stage %in% c("all", "finalize")) {
        if (is.null(output))
            cli::cli_abort("{.arg output} is required for stage {.val {stage}}.")
        .finalize_stage(tmp_dir, output, ppm, grouper,
                        fdr_thresh, cor_thresh, min_count,
                        ion_mode, ...)
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
.finalize_stage <- function(tmp_dir, output, ppm, grouper,
                            fdr_thresh, cor_thresh, min_count,
                            ion_mode, ...) {
    fd      <- readRDS(file.path(tmp_dir, "fd.rds"))
    batches <- list.files(tmp_dir, pattern = "^features_batch_.*\\.parquet$",
                          full.names = TRUE)

    if (length(batches)) {
        cli::cli_progress_step("Combining {length(batches)} filter batch{?s}")
        cleaned <- data.table::rbindlist(
            lapply(batches, arrow::read_parquet),
            fill = TRUE
        )
    } else {
        cli::cli_progress_step("Filtering pseudo-spectra ({nrow(fd)} features)")
        cleaned <- clean_pseudo_spectra(
            targets        = rownames(fd),
            file_directory = tmp_dir,
            fd             = fd,
            grouper        = grouper,
            fdr_thresh     = fdr_thresh,
            cor_thresh     = cor_thresh,
            min_count      = min_count,
            report         = FALSE,
            mgf_dir        = NULL,
            progress       = FALSE,
            ...
        )
    }

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
.filter_stage <- function(tmp_dir, batch_idx, n_batches, grouper,
                          fdr_thresh, cor_thresh, min_count, ...) {
    fd <- readRDS(file.path(tmp_dir, "fd.rds"))
    cli::cli_progress_step(
        "Filtering pseudo-spectra (batch {batch_idx}/{n_batches})")
    cleaned <- clean_pseudo_spectra(
        targets        = rownames(fd),
        file_directory = tmp_dir,
        fd             = fd,
        grouper        = grouper,
        fdr_thresh     = fdr_thresh,
        cor_thresh     = cor_thresh,
        min_count      = min_count,
        batch_idx      = batch_idx,
        n_batches      = n_batches,
        report         = FALSE,
        mgf_dir        = NULL,
        progress       = TRUE,
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

#' Build pseudo-spectra MGF from mzML file(s) against a target list
#'
#' User-facing entry point for the targeted workflow: given a CSV of
#' `mz`/`rt` targets and one-or-many mzML files, produce an MGF of
#' pseudo-spectra. Same `stage`/`sample_idx` shape as [process_xcms()]
#' so a DAG can parallelize per-sample; interactive users leave
#' `stage = "all"`.
#'
#' # On-disk contract (identical to [process_xcms()])
#'
#' \itemize{
#'   \item `<tmp_dir>/manifest.csv` — `sample_idx,filename`.
#'   \item `<tmp_dir>/fd.rds` — feature-definition frame derived from
#'     the target list (`feature_id`, `mzmed`, `rtmed`).
#'   \item `<tmp_dir>/sample_<idx>_pool.parquet`,
#'         `<tmp_dir>/sample_<idx>_pointers.parquet`.
#' }
#'
#' @param input Single mzML path, a directory of mzMLs, or a character
#'   vector of mzML paths.
#' @param output MGF path. Required for `stage = "all"` or `"finalize"`.
#' @param targetlist Path to CSV with `mz` and `rt` columns (optional
#'   `feature`), or a data.frame with the same shape.
#' @param tmp_dir Working directory. `NULL` only for `stage = "all"`
#'   (creates + cleans a `tempfile`); required otherwise.
#' @param stage One of `"all"`, `"prepare"`, `"sample"`, `"filter"`,
#'   `"finalize"`. `"filter"` is the batched scale-out path;
#'   `"finalize"` auto-aggregates batch parquets if present.
#' @param sample_idx,sample_file Sample addressing for
#'   `stage = "sample"`; idx is 1-based into `manifest`.
#' @param batch_idx,n_batches For `stage = "filter"`: split `targets`
#'   into `n_batches` contiguous chunks and process chunk `batch_idx`
#'   (1-based). One `filter` DAG job per batch; `finalize` unions them.
#' @param ppm Fragment/precursor tolerance in ppm.
#' @param rt_width Total width (seconds) of the RT window around each
#'   target's `rt` when extracting scans.
#' @param grouper Fragment-grouping function; `NULL` uses
#'   `ppm_upgma(ppm = ppm)`.
#' @param fdr_thresh,cor_thresh,min_count Passed to
#'   [clean_pseudo_spectra()] (permissive defaults).
#' @param ion_mode `"POSITIVE"` or `"NEGATIVE"`.
#' @param ... Forwarded to [clean_pseudo_spectra()].
#' @return Invisibly the `tmp_dir` path.
#' @family api
#' @export
process_targetlist <- function(
    input,
    output      = NULL,
    targetlist  = NULL,
    tmp_dir     = NULL,
    stage       = c("all", "prepare", "sample", "filter", "finalize"),
    sample_idx  = NULL,
    sample_file = NULL,
    ppm         = 20,
    rt_width    = 15,
    grouper     = NULL,
    fdr_thresh  = 0.1,
    cor_thresh  = 0.1,
    min_count   = 10L,
    batch_idx   = NULL,
    n_batches   = 1L,
    ion_mode    = "POSITIVE",
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

    if (stage %in% c("all", "prepare")) {
        if (is.null(targetlist))
            cli::cli_abort("{.arg targetlist} is required for stage {.val {stage}}.")
        .prepare_targetlist_stage(input, targetlist, tmp_dir)
    }

    if (stage == "sample") {
        sid <- .resolve_sample_id(sample_idx, sample_file, tmp_dir)
        .sample_targetlist_stage(tmp_dir, sid, ppm, rt_width)
    }

    if (stage == "all") {
        manifest <- utils::read.csv(file.path(tmp_dir, "manifest.csv"))
        cli::cli_progress_bar("Building per-sample parquets",
                              total = nrow(manifest))
        for (sid in manifest$sample_idx) {
            suppressMessages(.sample_targetlist_stage(tmp_dir, sid, ppm, rt_width))
            cli::cli_progress_update()
        }
        cli::cli_progress_done()
    }

    if (stage == "filter") {
        if (is.null(batch_idx))
            cli::cli_abort("stage = {.val filter} requires {.arg batch_idx}.")
        .filter_stage(tmp_dir, batch_idx, n_batches, grouper,
                      fdr_thresh, cor_thresh, min_count, ...)
    }

    if (stage %in% c("all", "finalize")) {
        if (is.null(output))
            cli::cli_abort("{.arg output} is required for stage {.val {stage}}.")
        .finalize_stage(tmp_dir, output, ppm, grouper,
                        fdr_thresh, cor_thresh, min_count,
                        ion_mode, ...)
    }

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
    tl <- tl[!is.na(tl$mz) & !is.na(tl$rt), ]
    if (max(tl$rt) < 60) {
        cli::cli_alert_info(
            "RT looks like minutes (max = {round(max(tl$rt), 2)}); converting to seconds.")
        tl$rt <- tl$rt * 60
    }
    if (!"feature" %in% colnames(tl))
        tl$feature <- sprintf("FT%0*d", nchar(nrow(tl)), seq_len(nrow(tl)))
    tl$feature <- as.character(tl$feature)
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
        message("Sample ", sample_idx, ": no targets fell in the ", rt_width, "s window.")
        return(invisible(NULL))
    }

    max_diff <- max(fd$mzmed, na.rm = TRUE) * ppm * 1e-6
    Ms2PseudoRaw(
        single_sample  = scan_info,
        spectra_object = sps,
        output_dir     = tmp_dir,
        maxDiff        = max_diff
    )
}


# ---- carpediams (unified dispatch) -------------------------------------------

#' Build pseudo-spectra MGF - unified entry point
#'
#' Auto-dispatches to [process_mzml()], [process_targetlist()], or
#' [process_xcms()] based on the shape of `input` and whether a
#' `targetlist` is supplied. Same defaults, thresholds, and stage
#' machinery as the routed function.
#'
#' # Routing rules
#'
#' \itemize{
#'   \item An `XcmsExperiment`/`XCMSnExp` object, or a path to one
#'     (`.rds` / `.qs` / `.rda`) -> [process_xcms()].
#'     `targetlist` is ignored (with a warning) — the xcms object
#'     carries its own feature definitions.
#'   \item Any mzML input **with** a `targetlist` -> [process_targetlist()].
#'   \item A single mzML **without** a `targetlist` -> [process_mzml()].
#'   \item A directory or a character vector of mzMLs without a
#'     `targetlist` -> error (ambiguous — supply a targetlist or run
#'     [process_mzml()] per file).
#' }
#'
#' @param input Path or object; see routing rules.
#' @param output MGF path.
#' @param targetlist Optional CSV path or data.frame with `mz`, `rt`
#'   columns. Forces the targetlist route when combined with mzML input.
#' @param tmp_dir Working directory. `NULL` (default) creates and cleans
#'   a `tempfile()`. Pass an explicit path for split-stage runs.
#' @param ... Forwarded to the routed function. See that function's
#'   docs for stage machinery, thresholds, and route-specific params.
#' @return Whatever the routed function returns (typically an invisible
#'   `data.table` of pseudo-spectra or the `tmp_dir` path).
#' @family api
#' @export
carpediams <- function(input, output, targetlist = NULL, tmp_dir = NULL, ...) {
    route <- .detect_route(input, targetlist)
    cli::cli_h1("carpediams route: {.strong {route}}")

    if (is.null(tmp_dir)) {
        tmp_dir <- tempfile("carpediams_")
        dir.create(tmp_dir, recursive = TRUE)
        on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    } else {
        dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    }

    switch(route,
        xcms       = process_xcms(
            input   = input,
            output  = output,
            tmp_dir = tmp_dir,
            ...),
        targetlist = process_targetlist(
            input      = input,
            output     = output,
            targetlist = targetlist,
            tmp_dir    = tmp_dir,
            ...),
        mzml       = process_mzml(
            path    = input,
            tmp_dir = tmp_dir,
            output  = output,
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
