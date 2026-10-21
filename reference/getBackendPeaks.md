# Extract per-feature chromPeaks rows from an XCMSnExp

For each feature in `fd`, the rows of `all_peaks` listed in its
`peakidx`.

## Usage

``` r
getBackendPeaks(object, fd, all_peaks, sample_filter = NULL)
```

## Arguments

- object:

  Unused.

- fd:

  `featureDefinitions()` table with a `peakidx` list column.

- all_peaks:

  `chromPeaks()` of the object `fd` came from. Must be the full,
  unsubsetted table, because `peakidx` indexes its rows.

- sample_filter:

  Optional sample indices (into `all_peaks[, "sample"]`) to keep.

## Value

Named list (by feature id) of peak matrices; empty matrices for features
with no peaks.
