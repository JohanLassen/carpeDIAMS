#' One DeconvERROR placeholder block (no file I/O)
#'
#' `BEGIN IONS ... END IONS` with no peaks, titled
#' `<feature_id> DeconvERROR:<error_code>`. Used by [export_mgf()] so that
#' a feature with no real spectrum still gets exactly one MGF entry,
#' carrying the reason in its title rather than being silently absent.
#'
#' @keywords internal
#' @noRd
.deconv_error_block <- function(fid, mz_target, rt_target, error_code, ion_mode) {
  c(
    "BEGIN IONS",
    sprintf("TITLE=%s DeconvERROR:%s", fid, error_code),
    if (!is.na(mz_target)) sprintf("PEPMASS=%.4f", mz_target) else NULL,
    if (!is.na(rt_target)) sprintf("RTINSECONDS=%.2f", rt_target) else NULL,
    "CHARGE=1+",
    "MSLEVEL=2",
    sprintf("IONMODE=%s", ion_mode),
    "END IONS",
    ""
  )
}


#' Write a single-feature MGF block (concat-friendly)
#'
#' Writes one `BEGIN IONS ... END IONS` block (`TITLE`, `PEPMASS`,
#' optional `RTINSECONDS`, `CHARGE=1+`, peaks). Used by
#' [clean_pseudo_spectra()] when `mgf_dir` is set. Files can be joined
#' with `cat`.
#'
#' @param pseudo_spec Filtered per-group table with `fragment_mz` and
#'   `int`.
#' @param mz_target Numeric. Precursor m/z.
#' @param feature_id Used as `TITLE=`.
#' @param path Output `.mgf` path.
#' @param rt Optional retention time (seconds); skipped when `NA`.
#' @return Invisibly `path`, or `NULL` when `pseudo_spec` is empty.
#' @keywords internal
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
 

#' Write a list of spectra as one MGF and optional SIRIUS .ms files
#'
#' Standalone writer for spectra you assembled yourself; the pipeline
#' uses [export_mgf()]. Writes in chunks of `chunk_size` blocks, and one
#' `feature_<id>.ms` file per spectrum when `sirius_dir` is set.
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
#' Writes one MGF block per row of `fd`, in `fd` order. Features with
#' fragment rows get `TITLE`, `PEPMASS` (from `mzmed`), `RTINSECONDS`
#' (from `rtmed`), `CHARGE=1+`, `MSLEVEL=2`, `IONMODE` and the peaks
#' sorted by m/z. Features with an error row, or missing from `results`,
#' get an empty block titled `<feature_id> DeconvERROR:<reason>`
#' (`NoAssociatedFragments` when missing). Warns if no feature has
#' peaks.
#'
#' @param results Output of [clean_pseudo_spectra()]: `feature_id`, `mz`,
#'   `total_intensity`, optional `deconv_error`.
#' @param fd Feature table with `mzmed`/`rtmed` (or `mz`/`rt`); rownames or
#'   a `feature_id` column give the ids.
#' @param mgf_path Output path (directories are created).
#' @param ion_mode Written as `IONMODE=`.
#' @return Invisibly `mgf_path`.
#' @export
export_mgf <- function(results, fd, mgf_path, ion_mode = "POSITIVE") {
  # --- precursor m/z + RT lookup from fd ---
  fd_df <- as.data.frame(fd)
  if ("feature_id" %in% colnames(fd_df)) rownames(fd_df) <- fd_df$feature_id
  if (!"mzmed" %in% colnames(fd_df) && "mz" %in% colnames(fd_df)) fd_df$mzmed <- fd_df$mz
  if (!"rtmed" %in% colnames(fd_df) && "rt" %in% colnames(fd_df)) fd_df$rtmed <- fd_df$rt

  results_dt <- if (is.null(results)) data.table::data.table() else data.table::as.data.table(results)
  has_feature_col <- "feature_id" %in% colnames(results_dt)
  has_error_col   <- "deconv_error" %in% colnames(results_dt)

  dir.create(dirname(mgf_path), showWarnings = FALSE, recursive = TRUE)
  con <- file(mgf_path, open = "w")
  on.exit(close(con), add = TRUE)

  n_real <- 0L
  for (fid in rownames(fd_df)) {
    mz_target <- fd_df[fid, "mzmed"]
    rt_target <- fd_df[fid, "rtmed"]

    block <- if (has_feature_col) results_dt[feature_id == fid] else results_dt[0]
    is_placeholder <- nrow(block) == 0 ||
      (has_error_col && !is.na(block$deconv_error[1]))

    if (is_placeholder) {
      reason <- if (nrow(block) && has_error_col) block$deconv_error[1] else "NoAssociatedFragments"
      writeLines(.deconv_error_block(fid, mz_target, rt_target, reason, ion_mode), con)
      next
    }

    n_real <- n_real + 1L
    block <- block[order(mz)]                       # fragments ascending m/z
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

  if (n_real == 0L)
    warning("No spectra to export - every feature was written as a DeconvERROR placeholder; see MGF titles for reasons.")
  invisible(mgf_path)
}