
get_raw_spec <- function(feature_id, file_directory, fd = NULL, con = NULL) {
  own_con <- is.null(con)
  if (own_con) {
    con <- setup_connection(file_directory)
    on.exit(dbDisconnect(con, shutdown = TRUE))
  }
  dt <- batch_get_feature_data(feature_id, con)

  # If fd is provided, filter using the true precursor mz from featureDefinitions.
  # This corrects for parquets built with the old getBackendPeaks, where mz_target
  # in the pointer files may have been taken from the wrong (foverlaps-matched) peak.
  if (!is.null(fd)) {
    fd_df   <- as.data.frame(fd)
    true_mz <- fd_df[as.character(feature_id), "mzmed"]
    dt <- dt[is.na(fragment_mz) | fragment_mz <= true_mz + 2]
  }

  dt[, feature_id := feature_id]
  dt
}



# ===========================================================================
# 6. FEATURE INSPECTION: Get raw MS1+MS2 data for a single feature
#    - feature_id: integer index (1070), XCMS name ("FT1070"), or plain string
#    - parquet_dir: path to pre-computed parquets (fast); omit to compute live
#    - msobject: XCMSnExp — only needed when parquet_dir is NULL
#    - con: persistent DuckDB connection (reuse across calls for speed)
#    Output columns match batch_get_feature_data:
#      sample, feature, mz, fragment_mz, int, ms1_scan, ms2_scan, ms1_int
# ===========================================================================
get_feature_data <- function(
  feature_id,
  fd,
  msobject    = NULL,
  parquet_dir = NULL,
  con         = NULL,
  maxDiff     = 0.01
) {
  fid <- as.character(feature_id)

  # --- Fast path: pre-computed parquets ---
  if (!is.null(parquet_dir)) {
    own_con <- is.null(con)
    if (own_con) {
      con <- setup_connection(parquet_dir)
      on.exit(dbDisconnect(con, shutdown = TRUE))
    }
    return(batch_get_feature_data(fid, con))
  }

  # --- Live path: extract directly from msobject ---
  if (is.null(msobject)) stop("Provide either `parquet_dir` or `msobject`")

  fd_df <- as.data.frame(fd)
  if (!fid %in% rownames(fd_df)) stop(paste("feature_id not found in fd:", fid))

  peak_indices <- fd_df$peakidx[[fid]]
  mz_target    <- fd_df$mzmed[[fid]]

  all_peaks <- as.data.table(chromPeaks(msobject))
  all_peaks[, global_idx := .I]
  feat_peaks <- all_peaks[global_idx %in% peak_indices]

  if (nrow(feat_peaks) == 0) return(data.table())

  filenames  <- fileNames(msobject)
  sp_obj     <- spectra(msobject)
  spectra_df <- sp_obj@backend@spectraData@listData |> bind_cols()
  sp_by_file <- split(spectra_df, spectra_df$dataOrigin)

  results <- lapply(seq_len(nrow(feat_peaks)), function(i) {
    pk    <- feat_peaks[i]
    fname <- filenames[pk$sample]

    file_sp  <- sp_by_file[[fname]]
    if (is.null(file_sp)) return(NULL)

    in_range <- file_sp$rtime_adjusted >= pk$rtmin &
                file_sp$rtime_adjusted <= pk$rtmax
    window   <- file_sp[in_range, ]
    if (nrow(window) == 0) return(NULL)

    spec <- sp_obj |>
      filterDataOrigin(fname) |>
      filterAcquisitionNum(as.integer(window$acquisitionNum))
    if (length(spec) == 0) return(NULL)

    pks_list <- peaksData(spec)
    n_pks    <- vapply(pks_list, nrow, integer(1L))
    pks_mat  <- do.call(rbind, pks_list)

    all_ms <- data.table(
      frag_mz = pks_mat[, 1],
      int     = pks_mat[, 2],
      scan    = rep(acquisitionNum(spec), n_pks),
      mslevel = rep(msLevel(spec), n_pks)
    )

    ms1 <- all_ms[mslevel == 1L, .(ms1_int = sum(int)), by = .(ms1_scan = scan)]
    ms1[, ms2_scan := ms1_scan + 1L]

    ms2 <- all_ms[
      mslevel == 2L & frag_mz <= mz_target + 2,
      .(fragment_mz = frag_mz, int, ms2_scan = scan)
    ]

    # LEFT JOIN so every MS1 row is kept (fragment_mz/int = NA when no MS2 follows)
    joined <- merge(ms1, ms2, by = "ms2_scan", all.x = TRUE)
    joined[, `:=`(sample = pk$sample, feature = fid, mz = mz_target)]
    joined[, .(sample, feature, mz, fragment_mz, int, ms1_scan, ms2_scan, ms1_int)]
  })

  rbindlist(Filter(Negate(is.null), results))
}
