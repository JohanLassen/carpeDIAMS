# Build a pseudo-spectra MGF from one mzML file

Single-sample route with no prior feature list. Steps:

## Usage

``` r
process_mzml(
  path,
  tmp_dir,
  output = NULL,
  rt_range = c(30, 1200),
  grouper = NULL,
  group_statistics_fn = NULL,
  filter_fn = NULL
)
```

## Arguments

- path:

  mzML file path.

- tmp_dir:

  Directory for the intermediate parquets.

- output:

  MGF path. Required: the function stops before reading any data if it
  is `NULL`.

- rt_range:

  Retention-time window (seconds) to keep.

- grouper:

  A grouper closure (see
  [grouper-contract](https://johanlassen.github.io/carpeDIAMS/reference/grouper-contract.md));
  `NULL` uses `ppm_upgma(ppm = 30)`.

- group_statistics_fn:

  A statistics closure (see
  [stats-contract](https://johanlassen.github.io/carpeDIAMS/reference/stats-contract.md));
  `NULL` uses `correlation_statistics(min_count = 10)`.

- filter_fn:

  A filter closure (see
  [filter-contract](https://johanlassen.github.io/carpeDIAMS/reference/filter-contract.md));
  `NULL` uses
  `correlation_filtering(cor_thresh = 0.1, fdr_thresh = 0.1)`.

## Value

The
[`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md)
table. The MGF is written as a side effect.

## Details

1.  Read the file and keep spectra inside `rt_range`.

2.  If no spectrum is labelled MS2, assign MS levels by scan-index
    parity (bbCID alternates MS1/MS2). Aborts if the split is badly
    unbalanced.

3.  Detect chromatographic peaks with centWave (fixed parameters:
    `peakwidth = c(2, 30)`, `snthresh = 6`, `ppm = 12`,
    `prefilter = c(3, 1000)`, `mzdiff = 0.01`).

4.  Write pool/pointer parquets to `tmp_dir` and build one
    pseudo-spectrum per peak with
    [`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md).

5.  Write the MGF with
    [`export_mgf()`](https://johanlassen.github.io/carpeDIAMS/reference/export_mgf.md).

## See also

Other api:
[`carpediams()`](https://johanlassen.github.io/carpeDIAMS/reference/carpediams.md),
[`process_targetlist()`](https://johanlassen.github.io/carpeDIAMS/reference/process_targetlist.md),
[`process_xcms()`](https://johanlassen.github.io/carpeDIAMS/reference/process_xcms.md)
