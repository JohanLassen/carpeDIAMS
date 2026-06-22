  #' Write pseudo-spectra as a combined MGF + per-feature SIRIUS .ms files
  #'
  #' Streams `spectra_list` to one MGF (chunked to avoid the 2^31-1 byte
  #' writeLines limit), and optionally writes one SIRIUS `.ms` file per
  #' feature into `sirius_dir`.
  #'
  #' Each element of `spectra_list` must be a list with:
  #'   \itemize{
  #'     \item `peaks` — data.frame with columns `mz`, `total_intensity`
  #'     \item `feature_id` — character, used as MGF `TITLE=` and SIRIUS `>compound`
  #'     \item `mz_target` — numeric, precursor m/z
  #'     \item `rtmed` — numeric, retention time in seconds
  #'     \item `ESI` — character, "POSITIVE" / "NEGATIVE"
  #'   }
  #'
  #' @param spectra_list As described above.
  #' @param mgf_path Output `.mgf` path.
  #' @param sirius_dir Optional directory for per-feature `.ms` files. `NULL`
  #'   (default) skips SIRIUS output entirely.
  #' @param chunk_size Integer. MGF entries buffered before each `writeLines`.
  #' @return Invisibly `NULL`.
  #' @export
  export_pseudospectra <- function(spectra_list,
                                   mgf_path,
                                   sirius_dir = NULL,
                                   chunk_size = 1000L) {

    if (!is.null(sirius_dir)) {
      dir.create(sirius_dir, showWarnings = FALSE, recursive = TRUE)
    }

    mgf_con <- file(mgf_path, open = "w")
    on.exit(close(mgf_con), add = TRUE)

    n       <- length(spectra_list)
    buf     <- vector("character", min(chunk_size, n))
    buf_idx <- 0L

    for (i in seq_along(spectra_list)) {
      spec       <- spectra_list[[i]]
      peak_data  <- spec$peaks[order(spec$peaks$mz), ]
      peak_block <- paste(peak_data$mz, peak_data$total_intensity,
                          sep = " ", collapse = "\n")

      buf_idx       <- buf_idx + 1L
      buf[buf_idx]  <- paste0(
        "BEGIN IONS\n",
        "TITLE=", spec$feature_id, "\n",
        "PEPMASS=", spec$mz_target, "\n",
        "CHARGE=1+\n",
        "MSLEVEL=2\n",
        "RTINSECONDS=", spec$rtmed, "\n",
        "IONMODE=", spec$ESI, "\n",
        peak_block, "\n",
        "END IONS\n"
      )

      if (buf_idx >= chunk_size) {
        writeLines(paste(buf[seq_len(buf_idx)], collapse = "\n"),
                   mgf_con, sep = "\n")
        buf_idx <- 0L
      }

      if (!is.null(sirius_dir)) {
        ms_content <- paste0(
          ">compound ", spec$feature_id, "\n",
          ">parentmz ", spec$mz_target, "\n",
          ">rt ",       spec$rtmed,     "\n",
          ">charge 1+\n",
          "\n>ms2\n",
          peak_block, "\n"
        )
        writeLines(
          ms_content,
          file.path(sirius_dir, paste0("feature_", spec$feature_id, ".ms"))
        )
      }
    }

    if (buf_idx > 0L) {
      writeLines(paste(buf[seq_len(buf_idx)], collapse = "\n"),
                 mgf_con, sep = "\n")
    }
    invisible(NULL)
  }