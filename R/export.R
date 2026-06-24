#' Write a single-feature MGF block (concat-friendly)
#'
#' Each file is a complete `BEGIN IONS ... END IONS` — merge with
#' `cat <dir>/*.mgf > combined.mgf`.
#'
#' @param pseudo_spec Output of [select_pseudo_spec()] (or a subset).
#' @param mz_target Numeric. Precursor m/z.
#' @param feature_id Used as `TITLE=`.
#' @param path Output `.mgf` path.
#' @param rt Optional retention time (seconds); skipped when `NA`.
#' @return Invisibly `path`, or `NULL` when `pseudo_spec` is empty.
#' @export
write_simple_mgf <- function(pseudo_spec, mz_target, feature_id, path,
                             rt = NA_real_) {
  if (is.null(pseudo_spec) || nrow(pseudo_spec) == 0) return(invisible(NULL))
  ord <- order(pseudo_spec$fragment_mz)
  body <- sprintf("%.4f %.0f", pseudo_spec$fragment_mz[ord], pseudo_spec$int[ord])
  header <- c(
    "BEGIN IONS",
    sprintf("TITLE=%s", feature_id),
    sprintf("PEPMASS=%.4f", mz_target),
    if (!is.na(rt)) sprintf("RTINSECONDS=%.2f", rt) else NULL,
    "CHARGE=1+"
  )
  writeLines(c(header, body, "END IONS", ""), path)
  invisible(path)
}
 

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


#' Assemble a combined MGF from a flat pseudo-spectra results table
#'
#' @param results data.table from clean_pseudo_spectra: one row per fragment,
#'   with columns feature_id, mz (fragment m/z), total_intensity.
#' @param fd featureDefinitions (or data.frame) giving precursor m/z + RT per
#'   feature. Rownames or a feature_id column must match results$feature_id.
#' @param mgf_path output .mgf path.
#' @param ion_mode "POSITIVE" or "NEGATIVE".
#' @return invisibly mgf_path.
#' @export
export_mgf <- function(results, fd, mgf_path, ion_mode = "POSITIVE") {
  if (is.null(results) || nrow(results) == 0) {
    warning("No spectra to export."); return(invisible(NULL))
  }

  # --- precursor m/z + RT lookup from fd ---
  fd_df <- as.data.frame(fd)
  if ("feature_id" %in% colnames(fd_df)) rownames(fd_df) <- fd_df$feature_id
  if (!"mzmed" %in% colnames(fd_df) && "mz" %in% colnames(fd_df)) fd_df$mzmed <- fd_df$mz
  if (!"rtmed" %in% colnames(fd_df) && "rt" %in% colnames(fd_df)) fd_df$rtmed <- fd_df$rt

  feats <- unique(results$feature_id)
  con <- file(mgf_path, open = "w")
  on.exit(close(con), add = TRUE)

  for (fid in feats) {
    block <- results[feature_id == fid]
    block <- block[order(mz)]                       # fragments ascending m/z

    has_fd    <- fid %in% rownames(fd_df)
    mz_target <- if (has_fd) fd_df[fid, "mzmed"] else NA_real_
    rt_target <- if (has_fd) fd_df[fid, "rtmed"] else NA_real_

    peaks <- sprintf("%.4f %.0f", block$mz, block$total_intensity)

    header <- c(
      "BEGIN IONS",
      sprintf("TITLE=%s", fid),
      if (!is.na(mz_target)) sprintf("PEPMASS=%.4f", mz_target) else NULL,
      if (!is.na(rt_target)) sprintf("RTINSECONDS=%.2f", rt_target) else NULL,
      "CHARGE=1+",
      "MSLEVEL=2",
      sprintf("IONMODE=%s", ion_mode)
    )
    writeLines(c(header, peaks, "END IONS", ""), con)
  }
  invisible(mgf_path)
}