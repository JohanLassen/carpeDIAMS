# Build a pseudo-spectra MGF from a grouped xcms object

Multi-sample route. Every feature in `featureDefinitions()` gets one MGF
block, built from its top-`top_n` most intense samples. Run it in one
call (`stage = "all"`) or one stage at a time so a scheduler can
parallelize the per-sample and per-batch work.

## Usage

``` r
process_xcms(
  input,
  output = NULL,
  tmp_dir = NULL,
  stage = c("all", "prepare", "sample", "filter", "finalize"),
  sample_idx = NULL,
  sample_file = NULL,
  subsample = FALSE,
  subset_size = 100L,
  q_prob = 0.9,
  top_n = 10L,
  ppm = 20,
  rt_buffer = 5,
  grouper = NULL,
  group_statistics_fn = NULL,
  filter_fn = NULL,
  fdr_thresh = 0.1,
  cor_thresh = 0.1,
  min_count = 10L,
  batch_idx = NULL,
  n_batches = 1L,
  ion_mode = "POSITIVE",
  ...
)
```

## Arguments

- input:

  An `XcmsExperiment`/`XCMSnExp` (after
  [`xcms::groupChromPeaks()`](https://rdrr.io/pkg/xcms/man/groupChromPeaks.html)),
  or a `.rds` / `.qs` / `.rda` path to one. Only read by `prepare`.

- output:

  MGF path. Required for `"all"` and `"finalize"`.

- tmp_dir:

  Working directory. `NULL` is only allowed for `stage = "all"`, which
  then uses a [`tempfile()`](https://rdrr.io/r/base/tempfile.html)
  deleted on exit.

- stage:

  One of `"all"`, `"prepare"`, `"sample"`, `"filter"`, `"finalize"`.

- sample_idx:

  Row of `manifest.csv` to process at `"sample"`. Takes precedence over
  `sample_file`.

- sample_file:

  mzML basename to process at `"sample"`; looked up in `manifest.csv`.

- subsample:

  If `TRUE`, keep `subset_size` samples chosen by apricot
  `MaxCoverageSelection` (needs Python with `apricot-select` and
  `numpy`). See
  [`subset_samples_xcms()`](https://johanlassen.github.io/carpeDIAMS/reference/subset_samples_xcms.md).

- subset_size, q_prob:

  Passed to
  [`subset_samples_xcms()`](https://johanlassen.github.io/carpeDIAMS/reference/subset_samples_xcms.md).

- top_n:

  Number of most intense samples per feature to extract. Used by
  `prepare` only.

- ppm:

  Used twice: at `sample`, the MS1 m/z tolerance for the precursor trace
  (converted to Da at the sample's largest feature m/z); and for the
  default grouper `ppm_upgma(ppm = ppm)`.

- rt_buffer:

  Seconds added to both sides of each peak's `rtmin`/`rtmax` at
  `sample`.

- grouper:

  A grouper closure (see
  [grouper-contract](https://johanlassen.github.io/carpeDIAMS/reference/grouper-contract.md));
  `NULL` uses `ppm_upgma(ppm = ppm)`.

- group_statistics_fn:

  A statistics closure (see
  [stats-contract](https://johanlassen.github.io/carpeDIAMS/reference/stats-contract.md));
  `NULL` uses `correlation_statistics(min_count = min_count)`.

- filter_fn:

  A filter closure (see
  [filter-contract](https://johanlassen.github.io/carpeDIAMS/reference/filter-contract.md));
  `NULL` uses `correlation_filtering(cor_thresh, fdr_thresh)`.

- fdr_thresh, cor_thresh, min_count:

  Only used to build the default `filter_fn` / `group_statistics_fn`.

- batch_idx, n_batches:

  For `"filter"`: process slice `batch_idx` (1-based) of `n_batches`. A
  slice past the last feature writes an empty parquet.

- ion_mode:

  `"POSITIVE"` or `"NEGATIVE"`, written as `IONMODE=`.

- ...:

  Forwarded to
  [`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md)
  at `"filter"` and `"finalize"` (e.g. `buffer`, `consensus_mz_fn`,
  `apex_source`, `window_mode`).

## Value

Invisibly `tmp_dir`.

## Stages

- `prepare`:

  Loads `input`, builds the feature × sample intensity matrix
  (`featureValues(method = "medret", value = "into")`), optionally
  subsamples with apricot, and ranks the top-`top_n` samples per feature
  (ties can add more). Writes `manifest.csv`, `mssubset.rds`, `fd.rds`,
  `sample_topfeatures.rds`.

- `sample`:

  For one sample: extracts the MS1/MS2 scans inside each qualifying
  feature's chromatographic peak (± `rt_buffer`) and writes
  `sample_<idx>_pool.parquet` (MS2) and `sample_<idx>_pointers.parquet`
  (MS1 trace per feature).

- `filter`:

  For one contiguous slice (`batch_idx` of `n_batches`) of the features:
  groups, scores and filters fragments with
  [`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md)
  and writes `features_batch_<batch_idx>.parquet`.

- `finalize`:

  If any `features_batch_*.parquet` exist in `tmp_dir`, binds **all** of
  them; otherwise runs the filtering for every feature in-process. Then
  writes the MGF with
  [`export_mgf()`](https://johanlassen.github.io/carpeDIAMS/reference/export_mgf.md).

- `all`:

  `prepare`, every `sample`, then `finalize` (in-process filtering, no
  batches).

All stages of one run must share `tmp_dir`.

## See also

Other api:
[`carpediams()`](https://johanlassen.github.io/carpeDIAMS/reference/carpediams.md),
[`process_mzml()`](https://johanlassen.github.io/carpeDIAMS/reference/process_mzml.md),
[`process_targetlist()`](https://johanlassen.github.io/carpeDIAMS/reference/process_targetlist.md)
