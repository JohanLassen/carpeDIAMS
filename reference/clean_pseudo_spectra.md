# Clean pseudo-spectra for many features (main batch entry)

Runs
[`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md)
and `filter_fn` for every target, reading the per-sample parquets in
`file_directory` through DuckDB. This is the work behind the `filter`
and `finalize` stages.

## Usage

``` r
clean_pseudo_spectra(
  targets,
  file_directory,
  fd = NULL,
  con = NULL,
  grouper = ppm_upgma(ppm = 30),
  buffer = 50L,
  fdr_thresh = 0.01,
  cor_thresh = 0.3,
  min_count = 10L,
  group_statistics_fn = NULL,
  filter_fn = NULL,
  consensus_mz_fn = NULL,
  apex_source = c("recompute", "upstream"),
  window_mode = c("buffer", "fwhm"),
  require_rt = FALSE,
  batch_idx = NULL,
  n_batches = 1L,
  report = TRUE,
  report_path = "pseudo_specs_report.pdf",
  mgf_dir = NULL,
  fetch_chunk_size = 50L,
  progress = TRUE
)
```

## Arguments

- targets:

  Character vector of feature ids.

- file_directory:

  Directory with the per-sample parquets.

- fd:

  Optional table with `mzmed` and `rtmed` (or `mz`/`rt`), with rownames
  or a `feature_id` column matching `targets`. Without it the precursor
  m/z comes from the parquets and RT is unknown.

- con:

  Optional open connection from
  [`setup_connection()`](https://johanlassen.github.io/carpeDIAMS/reference/setup_connection.md);
  it is left open.

- grouper:

  A grouper closure
  ([grouper-contract](https://johanlassen.github.io/carpeDIAMS/reference/grouper-contract.md)).

- buffer:

  Passed to
  [`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md).

- fdr_thresh, cor_thresh:

  Used only to build the default `filter_fn`.

- min_count:

  Used only to build the default `group_statistics_fn`.

- group_statistics_fn:

  A statistics closure
  ([stats-contract](https://johanlassen.github.io/carpeDIAMS/reference/stats-contract.md));
  `NULL` uses `correlation_statistics(min_count)`.

- filter_fn:

  A filter closure
  ([filter-contract](https://johanlassen.github.io/carpeDIAMS/reference/filter-contract.md));
  `NULL` uses `correlation_filtering(cor_thresh, fdr_thresh)`.

- consensus_mz_fn, apex_source, window_mode:

  Passed to
  [`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md).

- require_rt:

  If `TRUE`, a target without RT in `fd` is a `"MissingTargetInfo"`
  error. The target-list route sets this.

- batch_idx, n_batches:

  Process only slice `batch_idx` (1-based) of `targets` split into
  `n_batches` contiguous chunks. A slice past the end returns an empty
  table.

- report, report_path:

  Write diagnostic plots per feature to a PDF (needs ggplot2). **On by
  default**; the `process_*()` functions turn it off.

- mgf_dir:

  If set, also write one MGF file per feature here.

- fetch_chunk_size:

  Features per DuckDB query.

- progress:

  Show a progress bar.

## Value

`data.table` with one row per kept fragment: `mz_g2` and `mz` (both the
consensus fragment m/z), `grouped_mz` (group id), `total_intensity`
(mean per-cell intensity, not a sum), `n_samples`, `count`, `r`,
`t_stat`, `p_value`, `fdr`, `feature_id`, `deconv_error` (`NA` on real
rows). Column names are what
[`export_mgf()`](https://johanlassen.github.io/carpeDIAMS/reference/export_mgf.md)
reads.

## Details

Every target gives either fragment rows or exactly one error row with
`deconv_error` set to:

- `"MissingTargetInfo"`: no precursor m/z (or no RT when
  `require_rt = TRUE`).

- `"MissingMs1Signal"`: no pointer rows for the target, or no pointer
  parquets at all.

- `"NoAssociatedFragments"`: no group survived statistics or filtering,
  or
  [`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md)
  raised an error (errors are not shown).
