#' Extract per-feature chromPeaks rows from an XCMSnExp
#'
#' For each feature in `fd`, the rows of `all_peaks` listed in its
#' `peakidx`.
#'
#' @param object Unused.
#' @param fd `featureDefinitions()` table with a `peakidx` list column.
#' @param all_peaks `chromPeaks()` of the object `fd` came from. Must be
#'   the full, unsubsetted table, because `peakidx` indexes its rows.
#' @param sample_filter Optional sample indices (into `all_peaks[, "sample"]`)
#'   to keep.
#' @return Named list (by feature id) of peak matrices; empty matrices
#'   for features with no peaks.
#' @keywords internal
#' @export
getBackendPeaks <- function(object, fd, all_peaks, sample_filter = NULL) {
  # all_peaks must be the chromPeaks from the FULL (pre-subset) msobject that
  # fd was built from — subsetting the msobject changes row indices and breaks
  
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
  file_spectra <- spectra_list[[dataOrigin]]
  if (is.null(file_spectra)) return(NULL)

  rt_col   <- if ("rtime_adjusted" %in% colnames(file_spectra)) "rtime_adjusted" else "rtime"
  in_range <- file_spectra[[rt_col]] >= rtmin & file_spectra[[rt_col]] <= rtmax
  window   <- file_spectra[in_range, ]
  if (nrow(window) == 0) return(NULL)

  apex_idx <- which.min(abs(window[[rt_col]] - rt))

  list(
    sample          = sample,
    filename        = dataOrigin,
    mz              = mz,
    feature         = feature_id,
    scan_indices    = window$acquisitionNum,
    apex_scan_index = window$acquisitionNum[apex_idx],
    n_scans         = nrow(window)
  )
}

#' Map feature peaks to scan indices in the Spectra backend
#'
#' For each (feature, sample), takes the most intense peak (largest
#' `into`) and lists the scans whose retention time lies in
#' `[rtmin, rtmax]`, plus the scan closest to the peak `rt`. Uses
#' `rtime_adjusted` when present, else `rtime`.
#'
#' @param feature_peaks Output of [getBackendPeaks()].
#' @param filenames File paths indexed by the peaks' `sample` column,
#'   matching `spectra_df$dataOrigin`.
#' @param spectra_df Per-scan metadata with `dataOrigin`,
#'   `acquisitionNum` and `rtime`, e.g.
#'   `as.data.frame(Spectra::spectraData(sps))`.
#' @param target_list If `TRUE`, keep every peak per (feature, sample)
#'   instead of only the most intense.
#' @return Tibble with `sample, filename, mz, feature, scan_indices`
#'   (list), `apex_scan_index, n_scans`.
#' @keywords internal
#' @export
extractScanInfo <- function(feature_peaks, filenames, spectra_df, target_list = FALSE){

  spectra_list <- split(spectra_df, spectra_df$dataOrigin)

  if (target_list){
    out <- feature_peaks |>
      purrr::map_dfr(tibble::as_tibble, .id = "feature_id") |>
      dplyr::mutate(dataOrigin = filenames[sample]) |>
      dplyr::select(feature_id, dataOrigin, rtmin, rtmax, rt, sample, mz) |>
      purrr::pmap_dfr(get_scans_optimized, spectra_list, .progress = TRUE)
    return(out)
  }

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
#' Reads the scans of one file between the lowest and highest scan in
#' `single_sample` (± `buffer`) and writes two parquets to `output_dir`:
#'
#' * `sample_<id>_pool.parquet`: every MS2 centroid in that scan range
#'   (`sample, scan, rt, mz, int`).
#' * `sample_<id>_pointers.parquet`: for each feature and each MS1 scan
#'   within its scan range ± `buffer`, the summed intensity of MS1
#'   centroids within `maxDiff` Da of the precursor (`feature, scan,
#'   mz_target, apex_scan, ms1_int, mz, rt, sample`).
#'
#' @param single_sample One file's rows from [extractScanInfo()] (or the
#'   same columns built from a target list).
#' @param spectra_object `Spectra` object containing that file.
#' @param output_dir Output directory.
#' @param maxDiff MS1 m/z tolerance in Da.
#' @param buffer Extra scans on each side.
#' @param return_spectra If `TRUE`, also return the two tables.
#' @return `NULL`, or `list(pool, pointer)` when `return_spectra = TRUE`.
#'   Returns `NULL` without writing if no scans are in range.
#' @keywords internal
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
  pks_list <- Spectra::peaksData(b) |> as.list()
  n_peaks  <- vapply(pks_list, nrow, integer(1))
  pks_matrix <- do.call(rbind, pks_list)
  all_ms <- tibble::tibble(
    mz = pks_matrix[, 1],
    int = pks_matrix[, 2],
    scan = rep(Spectra::acquisitionNum(b), n_peaks),
    rt = rep(Spectra::rtime(b), n_peaks),
    mslevel = rep(Spectra::msLevel(b), n_peaks)
  )
  rm(b, pks_list, pks_matrix); gc()

  # 2. CREATE THE POOL: All unique MS2 fragments for this sample
  # We only save this ONCE per sample. No duplication.
  ms2_pool <- all_ms |>
    dplyr::filter(mslevel == 2) |>
    dplyr::mutate(sample = sample_id) |>
    dplyr::select(sample, scan, rt, mz, int) |>
    dplyr::arrange(scan)

  arrow::write_parquet(ms2_pool, file.path(output_dir, paste0("sample_", sample_id, "_pool.parquet")))

  # Keep MS1 in RAM for the loop
  ms1_all <- all_ms[all_ms$mslevel == 1, c("mz", "int", "scan", "rt")]
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

  # Summarize (rt is constant per scan - carried through, not aggregated)
  pointer_results <- matched[,
    .(ms1_int = sum(int), mz = mean(mz), rt = rt[1]),
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
#' The work behind `process_xcms(stage = "sample")`: looks up the
#' sample's peaks for `qualifying_features`, widens their RT windows by
#' `rt_buffer`, maps them to scans with [extractScanInfo()], and calls
#' [Ms2PseudoRaw()]. Writes nothing if the sample has no qualifying
#' features or no matching scans.
#'
#' @param mssubset Output of [subset_samples_xcms()].
#' @param fd `featureDefinitions()` table.
#' @param out_dir Output directory (created if missing).
#' @param sample_id Index into `fileNames(mssubset$msobject)`.
#' @param qualifying_features Feature ids to build for this sample.
#' @param ppm MS1 tolerance in ppm, converted to Da at the largest
#'   feature m/z in this sample. `NULL` uses 0.005 Da.
#' @param rt_buffer Seconds added to both sides of each peak's RT window.
#' @return `NULL`.
#' @keywords internal
#' @export
build_parquets_xcms <- function(
    mssubset,
    fd,
    out_dir,
    sample_id,
    qualifying_features,
    ppm       = NULL,
    rt_buffer = 0
) {

    if (!dir.exists(out_dir)) {
        message("Output directory ", out_dir, " does not exist. Creating...")
        dir.create(out_dir, recursive = TRUE)
    }

    if (length(qualifying_features) == 0) {
        message("Sample ", sample_id, ": no features to build after top_n filter.")
        return(invisible(NULL))
    }

    msobject <- mssubset$msobject
    sps_obj <- MsExperiment::spectra(msobject)
    spectra_df <- as.data.frame(Spectra::spectraData(sps_obj))
    filenames <- xcms::fileNames(msobject)

    fd_sub <- fd[qualifying_features, , drop = FALSE]
    full_idx <- match(basename(filenames[sample_id]), basename(mssubset$all_files))
    feature_peaks <- getBackendPeaks(msobject, all_peaks = mssubset$peaks, fd = fd_sub, sample_filter = full_idx)
    feature_peaks <- lapply(feature_peaks, function(m) {
      if (nrow(m)) m[, "sample"] <- sample_id
      m
      })
    if (rt_buffer > 0) {
        feature_peaks <- lapply(feature_peaks, function(m) {
            if (nrow(m) == 0) return(m)
            m[, "rtmin"] <- m[, "rtmin"] - rt_buffer
            m[, "rtmax"] <- m[, "rtmax"] + rt_buffer
            m
        })
    }
    scan_info <- extractScanInfo(feature_peaks, filenames, spectra_df)

    if (nrow(scan_info) == 0) {
        message("Sample ", sample_id, ": no features to build after top_n filter.")
        return(invisible(NULL))
    }

    current_spectra <- sps_obj |>
        Spectra::filterDataOrigin(scan_info$filename[1])

    rm(mssubset, msobject, sps_obj, spectra_df, filenames, feature_peaks, fd, fd_sub); gc()

    max_diff <- if (!is.null(ppm)) {
        max(scan_info$mz, na.rm = TRUE) * ppm * 1e-6
    } else 0.005

    Ms2PseudoRaw(
        single_sample  = scan_info,
        spectra_object = current_spectra,
        maxDiff        = max_diff,
        output_dir     = out_dir
    )
}