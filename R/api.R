

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

