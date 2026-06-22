  #' carpeDIAMS: DIA to pseudo-DDA spectra
  #'
  #' Turns DIA (bbCID) acquisitions into pseudo-DDA spectra by correlating MS1
  #' elution with MS2 fragment elution across many samples. Produces an MGF
  #' that downstream library/annotation tools consume.
  #'
  #' @keywords internal
  #' @useDynLib carpeDIAMS, .registration = TRUE
  #' @importFrom Rcpp sourceCpp
  "_PACKAGE"

utils::globalVariables(c(
  # data.table NSE
  "global_idx", "sample", "feature", "scan", "int", "mz", "mz_target",
  "mz_lo", "mz_hi", "scan_lo", "scan_hi", "apex_scan", "mz_start", "mz_end",
  "mslevel", "scan_indices", "apex_scan_index",
  # dplyr/tidyr NSE
  "feature_id", "name", "value", "dataOrigin", "rtmin", "rtmax", "rt", "into",
   # spectral_filtering
  "apex", "group", "fragment_mz", "ms1_int", "ms1_scan", "ms2_scan",
  "fdr", "p_value", "correlation", "count", "n_samples"
))