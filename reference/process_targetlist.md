# Build a pseudo-spectra MGF from mzML file(s) and a target list

Targeted route: no peak detection. Every target is extracted from every
file, in a fixed RT window around its `rt`. Same stages and `tmp_dir`
layout as
[`process_xcms()`](https://johanlassen.github.io/carpeDIAMS/reference/process_xcms.md),
except that `prepare` writes only `manifest.csv` and `fd.rds` (no
subsampling, no top-N ranking).

## Usage

``` r
process_targetlist(
  input,
  output = NULL,
  targetlist = NULL,
  tmp_dir = NULL,
  stage = c("all", "prepare", "sample", "filter", "finalize"),
  sample_idx = NULL,
  sample_file = NULL,
  ppm = 20,
  rt_width = 15,
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

  One mzML path, a directory of mzMLs (non-recursive), or a character
  vector of mzML paths. Only read by `prepare`.

- output:

  MGF path. Required for `"all"` and `"finalize"`.

- targetlist:

  CSV path or data.frame (see Details). Required for `"all"` and
  `"prepare"`.

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

- ppm:

  Used twice: at `sample`, the MS1 m/z tolerance for the precursor trace
  (converted to Da at the largest target m/z); and for the default
  grouper `ppm_upgma(ppm = ppm)`.

- rt_width:

  Total width (seconds) of the extraction window, centred on each
  target's `rt`.

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
  at `"filter"` and `"finalize"`. `require_rt` is always `TRUE` on this
  route.

## Value

Invisibly `tmp_dir`.

## Details

Target-list handling:

- Needs `mz` and `rt` columns; `feature` is optional (ids `FT01`,
  `FT02`, … are generated when missing).

- If the largest `rt` is below 60 it is taken as minutes and converted
  to seconds.

- Rows with missing `mz` or `rt` are kept and reported as
  `DeconvERROR:MissingTargetInfo` in the MGF.

At `sample`, files with no MS2 labels get MS levels by scan-index parity
(as in
[`process_mzml()`](https://johanlassen.github.io/carpeDIAMS/reference/process_mzml.md)).
Samples where no target falls within the window write no parquets;
`stage = "all"` reports them.

## See also

Other api:
[`carpediams()`](https://johanlassen.github.io/carpeDIAMS/reference/carpediams.md),
[`process_mzml()`](https://johanlassen.github.io/carpeDIAMS/reference/process_mzml.md),
[`process_xcms()`](https://johanlassen.github.io/carpeDIAMS/reference/process_xcms.md)
