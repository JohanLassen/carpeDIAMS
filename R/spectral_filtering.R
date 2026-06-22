#' Open a DuckDB connection with views over pool/pointer parquets
#'
#' Creates two SQL views (`pointers`, `pool`) on a directory of per-sample
#' parquets produced by [Ms2PseudoRaw()]. Caller is responsible for
#' `DBI::dbDisconnect(con, shutdown = TRUE)` — or pass the connection to
#' [clean_pseudo_spectra()] which closes it on exit.
#'
#' @param file_directory Directory containing `*_pool.parquet` + `*_pointers.parquet`.
#' @return A live DuckDB connection.
#' @export
setup_connection <- function(file_directory) {
  # Resolve to absolute path — DuckDB may not share R's working directory
  abs_dir <- normalizePath(file_directory, mustWork = TRUE)
  # Use a local temp dir for DuckDB spill files to avoid /tmp exhaustion on shared nodes
  tmp_dir <- file.path(abs_dir, ".duckdb_tmp")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  con <- DBI::dbConnect(duckdb(), config = list(temp_directory = tmp_dir))
  DBI::dbExecute(con, sprintf(
    "CREATE OR REPLACE VIEW pointers AS SELECT * FROM read_parquet('%s/*_pointers.parquet')", 
    abs_dir
  ))
  DBI::dbExecute(con, sprintf(
    "CREATE OR REPLACE VIEW pool AS SELECT * FROM read_parquet('%s/*_pool.parquet')", 
    abs_dir
  ))
  con
}

#' Fetch raw MS2 + MS1 rows for a vector of features from DuckDB
#'
#' Joins `pointers ⨝ pool` on `(sample, scan + 1)` with the precursor cut
#' `m.mz <= p.mz_target + 2`. Chunks the feature list to bound `IN (…)`.
#'
#' @param feature_ids Character vector of feature ids.
#' @param con A connection from [setup_connection()].
#' @param chunk_size Integer. Features per `IN`-clause.
#' @return `data.table` with columns
#'   `sample, feature, mz, fragment_mz, int, ms1_scan, ms2_scan, ms1_int`.
#' @export
batch_get_feature_data <- function(feature_ids, con, chunk_size = 500) {
  # DuckDB is fast at IN clauses but we chunk to avoid overly long SQL
  chunks <- split(feature_ids, ceiling(seq_along(feature_ids) / chunk_size))
  
  results <- data.table::rbindlist(lapply(chunks, function(ids) {
    id_list <- paste0("'", ids, "'", collapse = ", ")
    query <- sprintf("
      WITH feat_pointers AS (
        SELECT * FROM pointers WHERE feature IN (%s)
      )
      SELECT
        p.sample,
        p.feature,
        p.mz_target AS mz,
        m.mz AS fragment_mz,
        m.int,
        p.scan AS ms1_scan,
        m.scan AS ms2_scan,
        p.ms1_int
      FROM feat_pointers p
      LEFT JOIN pool m ON (
        p.scan + 1 = m.scan AND
        p.sample = m.sample
      )
      WHERE (m.mz <= p.mz_target + 2 OR m.mz IS NULL)
      ORDER BY p.feature, p.sample, p.scan
    ", id_list)
    data.table::as.data.table(dbGetQuery(con, query))
  }))
  
  return(results)
}



#' Build a feature's pseudo-spectrum from its raw MS1+MS2 data.table
#'
#' Per-sample apex window → two-pass ppm-UPGMA grouping → per-group
#' Pearson correlation of summed MS2 vs MS1 + BH FDR.
#'
#' @param dt Output of [batch_get_feature_data()] restricted to one feature.
#' @param mz_target Numeric. Precursor m/z (drives the ppm tolerance).
#' @param ms2_ppm Numeric. Fragment-mz grouping tolerance in ppm.
#' @param buffer Integer. Scan window around apex per sample.
#' @param min_count Integer. Minimum scan-count per group.
#' @return list(grouped_signal, correlations), or `NULL`.
#' @export
process_pseudo_spec <- function(
  dt,
  mz_target,
  grouper   = ppm_upgma(ppm = 30),
  buffer    = 50L,
  min_count = 10L
) {
  if (nrow(dt) == 0) return(NULL)
  dt <- dt[!is.na(int) & !is.na(fragment_mz)]
  if (nrow(dt) == 0) return(NULL)

  dt[, apex := ms1_scan[which.max(ms1_int)], by = sample]
  dt <- dt[abs(ms1_scan - apex) < buffer]
  if (nrow(dt) == 0) return(NULL)

  dt[, group := grouper(fragment_mz, target_mz = mz_target)]
  dt <- dt[, .(
    int         = sum(int),
    ms1_int     = mean(ms1_int),
    fragment_mz = mean(fragment_mz),
    mz          = mean(mz)
  ), by = .(group, sample, ms1_scan)]

  dt[, group := grouper(fragment_mz, target_mz = mz_target)]
  grouped_signal <- dt[, .(
    int         = sum(int),
    ms1_int     = mean(ms1_int),
    fragment_mz = mean(fragment_mz),
    mz          = mean(mz)
  ), by = .(group, sample, ms1_scan)]

  group_counts <- grouped_signal[, .(count = .N), by = group][count > min_count]
  if (nrow(group_counts) == 0) {
    return(list(grouped_signal = grouped_signal, correlations = data.table::data.table()))
  }

  cor_results <- grouped_signal[group %in% group_counts$group, {
    n <- .N
    if (n < 3L) {
      list(p_value = NA_real_, correlation = NA_real_)
    } else {
      xc <- int - mean(int)
      yc <- ms1_int - mean(ms1_int)
      denom <- sqrt(sum(xc * xc) * sum(yc * yc))
      if (denom == 0) {
        list(p_value = NA_real_, correlation = NA_real_)
      } else {
        r <- max(-1, min(1, sum(xc * yc) / denom))
        if (abs(r) >= 1) {
          list(p_value = 0, correlation = r)
        } else {
          tval <- r * sqrt((n - 2L) / (1 - r * r))
          list(p_value = 2 * stats::pt(-abs(tval), df = n - 2L), correlation = r)
        }
      }
    }
  }, by = group]

  group_summary <- grouped_signal[group %in% group_counts$group, .(
    fragment_mz = mean(fragment_mz),
    mz          = mean(mz),
    int         = mean(int),
    ms1_int     = mean(ms1_int)
  ), by = group]

  correlations <- merge(group_summary, cor_results, by = "group")
  correlations <- merge(correlations, group_counts, by = "group")
  correlations[, fdr := stats::p.adjust(p_value, method = "fdr")]

  list(grouped_signal = grouped_signal, correlations = correlations)
}

#' Apply self-peak filter and top-N selection to per-group correlations
#'
#' Removes the precursor self-peak (within `ms1_ppm`, skipped when `<= 0`),
#' then keeps the top `cor_count` groups by correlation, then the top
#' `fdr_count` of those by FDR — matches the notebook's two-stage ranking.
#'
#' @param correlations `data.table` from `process_pseudo_spec()$correlations`.
#' @param mz_target Numeric. Precursor m/z.
#' @param ms1_ppm Numeric. Tolerance for self-peak filter; `<= 0` skips.
#' @param cor_count Integer. Keep top-N by correlation.
#' @param fdr_count Integer. Then keep top-N of those by FDR.
#' @return A `data.table` (possibly empty).
#' @export
select_pseudo_spec <- function(correlations, mz_target, ms1_ppm, cor_count, fdr_count) {
  if (nrow(correlations) == 0) return(correlations)
  out <- correlations
  if (!is.na(ms1_ppm) && ms1_ppm > 0) {
    out <- out[abs(fragment_mz - mz_target) / mz_target * 1e6 > ms1_ppm]
  }
  if (nrow(out) == 0) return(out)
  out <- out[order(-correlation)][seq_len(min(cor_count, .N))]
  out <- out[order(fdr)][seq_len(min(fdr_count, .N))]
  out
}

#' Internal: diagnostic plots for one pseudo-spectrum (ggplot2-gated).
#' Returns `list()` when ggplot2 is unavailable.
#' @keywords internal
#' @noRd
.pseudo_spec_plots <- function(grouped_signal, correlations, pseudo_spec, feature_id) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(list())
  q95 <- if (nrow(correlations)) stats::quantile(correlations$int, 0.95, na.rm = TRUE) else 0

  p_corr <- ggplot2::ggplot(correlations,
                            ggplot2::aes(x = correlation, y = fdr,
                                         color = int > q95)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::labs(title = paste0(feature_id, ": correlation vs FDR"),
                  color = "high int") +
    ggplot2::theme_minimal(base_size = 9)

  p_spec <- ggplot2::ggplot(pseudo_spec) +
    ggplot2::geom_segment(ggplot2::aes(x = fragment_mz, xend = fragment_mz,
                                       y = 0, yend = int, color = 1 - fdr))
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_spec <- p_spec + ggrepel::geom_text_repel(
      data = pseudo_spec,
      ggplot2::aes(x = fragment_mz, y = int, label = round(fragment_mz, 3)),
      direction = "x", segment.size = 0.3, size = 2.5
    )
  }
  p_spec <- p_spec +
    ggplot2::labs(title = paste0(feature_id, ": pseudo-spectrum"),
                  x = "fragment m/z", y = "intensity") +
    ggplot2::theme_minimal(base_size = 9)

  filtered <- grouped_signal[group %in% pseudo_spec$group]

  p_ms1 <- if (nrow(filtered)) {
    ms1_summary <- filtered[, .(ms1_int = mean(ms1_int)), by = .(sample, ms1_scan)]
    ggplot2::ggplot(ms1_summary, ggplot2::aes(x = ms1_scan, y = ms1_int)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ sample, ncol = 2, scales = "free") +
      ggplot2::labs(title = paste0(feature_id, ": MS1 chromatograms")) +
      ggplot2::theme_minimal(base_size = 8)
  } else NULL

  p_eic <- if (nrow(filtered)) {
    eic <- data.table::copy(filtered)[, int := scale(int)[, 1], by = .(sample, group)]
    ggplot2::ggplot(eic, ggplot2::aes(x = ms1_scan, y = int,
                                      color = factor(group))) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ sample, ncol = 2, scales = "free") +
      ggplot2::labs(title = paste0(feature_id, ": fragment EICs (scaled)")) +
      ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(legend.position = "none")
  } else NULL

  list(p_corr = p_corr, p_spec = p_spec, p_ms1 = p_ms1, p_eic = p_eic)
}

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

#' Clean pseudo-spectra for many features (main batch entry)
#'
#' Orchestrates: connect → fetch → per-feature [process_pseudo_spec()] →
#' [select_pseudo_spec()] → optional plot/MGF. Returns a `data.table` in
#' the schema downstream `03_make_mgf.R` expects.
#'
#' @param targets Character vector of feature ids.
#' @param file_directory Path to pool/pointer parquets.
#' @param fd Optional `featureDefinitions()`-shaped frame; supplies
#'   `mzmed` and `rtmed` per target.
#' @param con Optional pre-opened connection (else opened + closed internally).
#' @param ms1_ppm,ms2_ppm Numeric ppm tolerances.
#' @param buffer,fdr_count,cor_count,min_count Tunables passed downstream.
#' @param report,report_path PDF diagnostic report on/off + path.
#' @param mgf_dir Optional directory for per-feature MGFs.
#' @param fetch_chunk_size Features per DuckDB query.
#' @param progress Print a per-feature progress counter.
#' @return A `data.table` matching the legacy 02_filter_features schema.
#' @export
clean_pseudo_spectra <- function(
  targets,
  file_directory,
  fd               = NULL,
  con              = NULL,
  ms1_ppm          = 0,
  grouper          = ppm_upgma(ppm = 30),
  buffer           = 50L,
  fdr_count        = 30L,
  cor_count        = 100L,
  min_count        = 10L,
  report           = TRUE,
  report_path      = "pseudo_specs_report.pdf",
  mgf_dir          = NULL,
  fetch_chunk_size = 200L,
  progress         = TRUE
) {
  own_con <- is.null(con)
  if (own_con) {
    con <- setup_connection(file_directory)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  if (report) {
    dir.create(dirname(report_path), showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(report_path, width = 11, height = 8)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  if (!is.null(mgf_dir)) dir.create(mgf_dir, showWarnings = FALSE, recursive = TRUE)

  fd_df <- NULL
  if (!is.null(fd)) {
    fd_df <- as.data.frame(fd)
    if ("feature_id" %in% colnames(fd_df)) rownames(fd_df) <- fd_df$feature_id
  }
  has_fd <- !is.null(fd_df)

  total <- length(targets)
  done  <- 0L
  fetch_chunks <- split(targets, ceiling(seq_along(targets) / fetch_chunk_size))

  results <- data.table::rbindlist(lapply(fetch_chunks, function(chunk_ids) {
    raw <- batch_get_feature_data(chunk_ids, con, chunk_size = length(chunk_ids))
    if (nrow(raw) == 0) {
      done <<- done + length(chunk_ids)
      return(NULL)
    }
    by_feature <- split(raw, by = "feature")

    data.table::rbindlist(lapply(chunk_ids, function(target) {
      done <<- done + 1L
      if (progress) cat(sprintf("\r  [%d/%d] %s         ", done, total, target))

      dt <- by_feature[[target]]
      if (is.null(dt) || nrow(dt) == 0) return(NULL)

      mz_target <- if (has_fd && target %in% rownames(fd_df)) fd_df[target, "mzmed"]
                   else dt$mz[1]
      if (length(mz_target) != 1 || is.na(mz_target)) return(NULL)
      rt_target <- if (has_fd && target %in% rownames(fd_df)) fd_df[target, "rtmed"]
                   else NA_real_

      out <- tryCatch(
        process_pseudo_spec(dt, mz_target = mz_target,
                            grouper = grouper, buffer = buffer,
                            min_count = min_count),
        error = function(e) NULL
      )
      if (is.null(out) || nrow(out$correlations) == 0) return(NULL)

      pseudo_spec <- select_pseudo_spec(out$correlations, mz_target,
                                        ms1_ppm, cor_count, fdr_count)
      if (nrow(pseudo_spec) == 0) return(NULL)

      if (report) {
        plots <- Filter(Negate(is.null),
                        .pseudo_spec_plots(out$grouped_signal, out$correlations,
                                           pseudo_spec, target))
        if (length(plots)) {
          if (requireNamespace("cowplot", quietly = TRUE)) {
            print(cowplot::plot_grid(plotlist = plots, ncol = 2))
          } else {
            for (p in plots) print(p)
          }
        }
      }

      if (!is.null(mgf_dir)) {
        write_simple_mgf(
          pseudo_spec, mz_target = mz_target, feature_id = target,
          path = file.path(mgf_dir, paste0(target, ".mgf")),
          rt = rt_target
        )
      }

      n_samples_per_group <- out$grouped_signal[
        group %in% pseudo_spec$group,
        .(n_samples = data.table::uniqueN(sample)),
        by = group
      ]
      pseudo_spec <- merge(pseudo_spec, n_samples_per_group,
                           by = "group", all.x = TRUE)

      data.table::data.table(
        mz_g2           = pseudo_spec$fragment_mz,
        grouped_mz      = as.integer(pseudo_spec$group),
        mz              = pseudo_spec$fragment_mz,
        total_intensity = pseudo_spec$int,
        n_samples       = as.integer(pseudo_spec$n_samples),
        count           = as.integer(pseudo_spec$count),
        r               = pseudo_spec$correlation,
        t_stat          = NA_real_,
        p_value         = pseudo_spec$p_value,
        fdr             = pseudo_spec$fdr,
        feature_id      = target
      )
    }), fill = TRUE)
  }), fill = TRUE)

  if (progress) cat("\n")
  results
}
