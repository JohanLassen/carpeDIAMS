#' Extract per-feature chromPeaks rows from an XCMSnExp
#'
#' Returns, per feature in `fd`, the matrix of `chromPeaks` rows that belong
#' to it (via `peakidx`). `all_peaks` must be the chromPeaks from the parent
#' (pre-subset) msobject — subsetting changes row indices and breaks the
#' `peakidx` <-> `chromPeaks` invariant. Pass `mssubset$peaks` (from
#' [subset_samples_xcms()]).
#'
#' @param object An `XCMSnExp` (only used for the default `all_peaks`).
#' @param fd `featureDefinitions()` table (rows = features, `peakidx` list-col).
#' @param all_peaks chromPeaks matrix from the parent msobject.
#' @param sample_filter Optional integer vector; keep only peaks from these samples.
#' @return Named list of chromPeaks-row matrices (one per feature id).
#' @export
getBackendPeaks <- function(object, fd, sample_filter = NULL) {
  # all_peaks must be the chromPeaks from the FULL (pre-subset) msobject that
  # fd was built from — subsetting the msobject changes row indices and breaks
  # the peakidx ↔ chromPeaks invariant. Pass mssubset$peaks here.
  
  all_peaks <- xcms::chromPeaks(object)
  peak_cols <- colnames(all_peaks)
  empty_mat <- matrix(nrow = 0, ncol = length(peak_cols),
                      dimnames = list(NULL, peak_cols))

  peaks_dt <- data.table::as.data.table(all_peaks)
  peaks_dt[, global_idx := .I]
  data.table::setkey(peaks_dt, global_idx)

  if (!is.null(sample_filter)) {
    peaks_dt <- peaks_dt[sample %in% sample_filter]
  }

  fd_df <- as.data.frame(fd)

  result <- lapply(fd_df$peakidx, function(idx) {
    if (is.null(idx) || length(idx) == 0) return(empty_mat)
    sub <- peaks_dt[global_idx %in% idx]
    if (nrow(sub) == 0) return(empty_mat)
    as.matrix(sub[, ..peak_cols])
  })
  names(result) <- rownames(fd_df)
  result
}

#' Internal: extract scan window for one (file, RT) tuple.
#' @keywords internal
#' @noRd
get_scans_optimized <- function(feature_id, dataOrigin, rtmin, rtmax, rt, sample, mz, spectra_list) {
  # Fast lookup in the pre-split list
  file_spectra <- spectra_list[[dataOrigin]]

  if (is.null(file_spectra)) return(NULL) # Safety check

  # Filter RT range (Vectorized)
  in_range <- file_spectra$rtime_adjusted >= rtmin &
    file_spectra$rtime_adjusted <= rtmax

  window <- file_spectra[in_range, ]

  if (nrow(window) == 0) return(NULL)

  # Find apex index efficiently
  apex_idx <- which.min(abs(window$rtime_adjusted - rt))

  return(
    list(
      sample = sample,
      filename = dataOrigin,
      mz = mz,
      # rt_apex_adjusted = rt,
      # rtmin = rtmin,
      # rtmax = rtmax,
      feature = feature_id,
      scan_indices = window$acquisitionNum,
      apex_scan_index = window$acquisitionNum[apex_idx],
      n_scans = nrow(window)
    )
  )
}

#' Map feature peaks to scan indices in the Spectra backend
#'
#' Joins per-feature chromPeaks (from [getBackendPeaks()]) to the Spectra
#' backend's per-scan metadata, producing one row per (feature, sample) with
#' the scan indices inside the peak RT window plus the apex scan.
#'
#' @param feature_peaks Output of [getBackendPeaks()].
#' @param filenames Character vector in `sample` order — typically
#'   `xcms::fileNames(msobject)`.
#' @param spectra_df Backend metadata, e.g.
#'   `Spectra::spectra(msobject)@backend@spectraData@listData |> dplyr::bind_cols()`.
#' @return Tibble: `sample, filename, mz, feature, scan_indices, apex_scan_index, n_scans`.
#' @export
extractScanInfo <- function(feature_peaks, filenames, spectra_df){

  spectra_list <- split(spectra_df, spectra_df$dataOrigin)

  out <- feature_peaks |>
     purrr::map_dfr(tibble::as_tibble, .id = "feature_id") |>
     dplyr::mutate(dataOrigin = filenames[sample]) |>
     dplyr::group_by(feature_id, sample) |>
     dplyr::slice_max(into, n = 1, with_ties = FALSE) |>
     dplyr::ungroup() |>
     dplyr::select(feature_id, dataOrigin, rtmin, rtmax, rt, sample, mz) |>
     purrr::pmap_dfr(get_scans_optimized, spectra_list, .progress = TRUE)

  return(out)
}

#' Build pool + pointer parquets for one sample
#'
#' Writes `sample_<id>_pool.parquet` (all MS2 fragments in the file) and
#' `sample_<id>_pointers.parquet` (one row per feature x MS1 scan, with
#' summed MS1 intensity). Together these feed [clean_pseudo_spectra()].
#'
#' @param single_sample Tibble from [extractScanInfo()] restricted to one file.
#' @param spectra_object Full `Spectra` object (`Spectra::spectra(msobject)`).
#' @param output_dir Directory to write the two parquets.
#' @param maxDiff Numeric. m/z tolerance for the MS1 ↔ precursor foverlaps join.
#' @param buffer Integer. Scan-window slack on each side of the feature's range.
#' @param return_spectra Logical. If `TRUE`, also return pool + pointers as tibbles.
#' @return `NULL` (or `list(pool, pointer)` if `return_spectra = TRUE`).
#' @export
Ms2PseudoRaw <- function(
    single_sample,
    spectra_object,
    output_dir = "results/",
    maxDiff = 0.01,
    buffer = 50L,
    return_spectra = FALSE
    ) {

  filename <- single_sample$filename[1]
  sample_id <- single_sample$sample[1]
  min_max_scan <- c(-buffer, buffer) + range(single_sample$scan_indices)
  # 1. Extraction
  b <- spectra_object |>
    Spectra::filterDataOrigin(filename) |>
    Spectra::filterAcquisitionNum(min_max_scan[1]:min_max_scan[2])
  rm(spectra_object); gc()

  if (length(b) == 0) return(NULL)

  pks_list <- Spectra::peaksData(b)
  n_peaks  <- vapply(pks_list, nrow, integer(1))
  pks_matrix <- do.call(rbind, pks_list)
  all_ms <- tibble::tibble(
    mz = pks_matrix[, 1],
    int = pks_matrix[, 2],
    scan = rep(Spectra::acquisitionNum(b), n_peaks),
    mslevel = rep(Spectra::msLevel(b), n_peaks)
  )
  rm(b, pks_list, pks_matrix); gc()

  # 2. CREATE THE POOL: All unique MS2 fragments for this sample
  # We only save this ONCE per sample. No duplication.
  ms2_pool <- all_ms |>
    dplyr::filter(mslevel == 2) |>
    dplyr::mutate(sample = sample_id) |>
    dplyr::select(sample, scan, mz, int) |>
    dplyr::arrange(scan)

  arrow::write_parquet(ms2_pool, file.path(output_dir, paste0("sample_", sample_id, "_pool.parquet")))

  # Keep MS1 in RAM for the loop
  ms1_all <- all_ms[all_ms$mslevel == 1, c("mz", "int", "scan")]
  rm(all_ms); gc()

  # 3. CREATE POINTERS: Loop through features
  feature_list <- split(single_sample, single_sample$feature)

  # Prepare feature ranges
  feature_ranges <- single_sample |>
    dplyr::group_by(feature) |>
    dplyr::summarize(
      mz_target = mz[1],
      mz_lo = mz[1] - maxDiff,
      mz_hi = mz[1] + maxDiff,
      scan_lo = min(scan_indices) - buffer,
      scan_hi = max(scan_indices) + buffer,
      apex_scan = apex_scan_index[1],
      .groups = "drop"
    ) |>
    data.table::as.data.table()

  feature_ranges[, `:=`(scan_lo_num = as.numeric(scan_lo), 
                        scan_hi_num = as.numeric(scan_hi))]
  data.table::setkey(feature_ranges, mz_lo, mz_hi)   
  
  ms1_dt <- data.table::as.data.table(ms1_all)
  ms1_dt[, `:=`(mz_start = mz, mz_end = mz)]
  # no setkey on ms1_dt!

  # First: overlap on mz interval
  matched <- data.table::foverlaps(ms1_dt, feature_ranges,
                      by.x = c("mz_start", "mz_end"),
                      by.y = c("mz_lo", "mz_hi"),
                      nomatch = NULL)

  # Then: filter scan range in a second step
  matched <- matched[scan >= scan_lo & scan <= scan_hi]

  # Summarize
  pointer_results <- matched[,
    .(ms1_int = sum(int), mz = mean(mz)),
    by = .(feature, scan, mz_target, apex_scan)
  ][, sample := sample_id] |> bind_rows()

  # Save the Pointers
  arrow::write_parquet(pointer_results, file.path(output_dir, paste0("sample_", sample_id, "_pointers.parquet")))

  rm(ms1_dt, matched); gc()
  message("Finished Sample ", sample_id)

  if (return_spectra) {
    return(list(
      pool = ms2_pool,
      pointer = pointer_results
    ))
  }
}


#' Build pool + pointer parquets for one sample, from a subsetted msobject
#'
#' Stage-01 orchestrator: takes the output of [subset_samples_xcms()] and
#' calls [Ms2PseudoRaw()] for one `sample_id`.
#'
#' @param mssubset Output of [subset_samples_xcms()].
#' @param fd `featureDefinitions()` table.
#' @param out_dir Directory for the two parquets.
#' @param sample_id Integer. 1-based index into the subsetted msobject.
#' @param include_quantile Currently unused; top-10 samples per feature
#'   are picked by `slice_max(n = 10)`.
#' @return Invisibly `NULL`. Writes two parquets to `out_dir`.
#' @export
build_parquets_xcms <- function(
    mssubset,
    fd,
    out_dir,
    sample_id,
    include_quantile = 0.90
) {
    
    # Check if output directory exists, if not create it
    if (!dir.exists(out_dir)) {
        message("Output directory ", out_dir, " does not exist. Creating...")
        dir.create(out_dir, recursive = TRUE)
    }

    # Load data from the subset
    msobject <- mssubset$msobject
    assay <- mssubset$assay
    spectra_df <- Spectra::spectra(msobject)@backend@spectraData@listData |> bind_cols()
    filenames <- xcms::fileNames(msobject)
    sample_file <- basename(filenames[sample_id]) 
    assay <- tibble::as_tibble(assay, rownames = "name")
    colnames(assay) <- c("name", rownames(fd))

    topX <- 
        assay |> 
        tidyr::pivot_longer(cols = -"name", names_to = "feature_id") |> 
        tidyr::replace_na(list(value = 0)) |> 
        dplyr::group_by(feature_id) |> 
        #dplyr::filter(value > quantile(value, include_quantile)) |> 
        dplyr::slice_max(value, n=10) |> 
        dplyr::ungroup()  |> 
        dplyr::filter(name == sample_file) 

    # Make scan info
    feature_peaks <- getBackendPeaks(msobject, fd = fd, all_peaks = mssubset$peaks, sample_filter = sample_id)
    scan_info <- extractScanInfo(feature_peaks, filenames, spectra_df) |> dplyr::filter(feature %in% topX$feature_id)
    
    current_spectra <- Spectra::spectra(msobject) |> Spectra::filterDataOrigin(scan_info$filename[1])
    
    # Clean up large objects
    rm(mssubset, msobject, assay, spectra_df, filenames, feature_peaks, fd, topX); gc()

    Ms2PseudoRaw(
        single_sample = scan_info,
        spectra_object = current_spectra,
        maxDiff = 0.005,
        output_dir = out_dir
    )
}