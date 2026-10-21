# Build a feature's per-group statistics table from raw MS1+MS2 data.table

Turns one feature's raw rows into a per-group statistics table:

## Usage

``` r
process_pseudo_spec(
  dt,
  mz_target,
  grouper = ppm_upgma(ppm = 30),
  group_statistics_fn = correlation_statistics(),
  buffer = 50L,
  consensus_mz_fn = NULL,
  apex_source = c("recompute", "upstream"),
  window_mode = c("buffer", "fwhm")
)
```

## Arguments

- dt:

  One feature's rows from
  [`batch_get_feature_data()`](https://johanlassen.github.io/carpeDIAMS/reference/batch_get_feature_data.md).

- mz_target:

  Precursor m/z, passed to the grouper as `target_mz`.

- grouper:

  A grouper closure
  ([grouper-contract](https://johanlassen.github.io/carpeDIAMS/reference/grouper-contract.md)).

- group_statistics_fn:

  A statistics closure
  ([stats-contract](https://johanlassen.github.io/carpeDIAMS/reference/stats-contract.md)).

- buffer:

  With `window_mode = "buffer"`, keep scans with
  `abs(ms1_scan - apex) < buffer`.

- consensus_mz_fn:

  A consensus-mz closure
  ([consensus-mz-contract](https://johanlassen.github.io/carpeDIAMS/reference/consensus-mz-contract.md));
  `NULL` uses
  [`consensus_mz_weighted()`](https://johanlassen.github.io/carpeDIAMS/reference/consensus_mz_weighted.md).

- apex_source:

  `"recompute"`: the scan with the highest MS1 intensity in `dt`.
  `"upstream"`: `dt$apex_scan`; samples where it is `NA` are dropped
  with a warning.

- window_mode:

  `"buffer"`: fixed scan window (see `buffer`). `"fwhm"`: walk out from
  the apex until MS1 intensity drops below half the apex intensity,
  within the scans present in `dt`.

## Value

The statistics table; zero rows if nothing is left at any step.

## Details

1.  Drop rows without a fragment.

2.  Find each sample's apex scan (`apex_source`) and keep the scans
    around it (`window_mode`).

3.  Group the fragment m/z of all samples together with `grouper`.

4.  Compute one consensus m/z per group (`consensus_mz_fn`).

5.  Score the groups with `group_statistics_fn`.

No filtering is applied;
[`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md)
adds that.
