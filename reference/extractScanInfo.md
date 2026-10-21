# Map feature peaks to scan indices in the Spectra backend

For each (feature, sample), takes the most intense peak (largest `into`)
and lists the scans whose retention time lies in `[rtmin, rtmax]`, plus
the scan closest to the peak `rt`. Uses `rtime_adjusted` when present,
else `rtime`.

## Usage

``` r
extractScanInfo(feature_peaks, filenames, spectra_df, target_list = FALSE)
```

## Arguments

- feature_peaks:

  Output of
  [`getBackendPeaks()`](https://johanlassen.github.io/carpeDIAMS/reference/getBackendPeaks.md).

- filenames:

  File paths indexed by the peaks' `sample` column, matching
  `spectra_df$dataOrigin`.

- spectra_df:

  Per-scan metadata with `dataOrigin`, `acquisitionNum` and `rtime`,
  e.g. `as.data.frame(Spectra::spectraData(sps))`.

- target_list:

  If `TRUE`, keep every peak per (feature, sample) instead of only the
  most intense.

## Value

Tibble with `sample, filename, mz, feature, scan_indices` (list),
`apex_scan_index, n_scans`.
