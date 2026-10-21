# Per-group correlation + apex-agreement statistics

As
[`correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_statistics.md),
plus two columns that compare where the fragment and the precursor
elute. Per sample, `m1` is the MS1-intensity-weighted mean scan of the
precursor trace, `m2` the fragment-intensity-weighted mean scan of the
group, and the shift is `m2 - m1` (in scans):

## Usage

``` r
apex_correlation_statistics(min_count = 10L)
```

## Arguments

- min_count:

  Groups need more than this many `(sample, ms1_scan)` cells to be
  scored.

## Value

A statistics closure — see
[stats-contract](https://johanlassen.github.io/carpeDIAMS/reference/stats-contract.md).

## Details

- `apex_shift_med`: median shift across samples.

- `apex_shift_iqr`: IQR of the shift across samples.

A real fragment co-elutes with its precursor, so its shift is near 0 in
every sample. Conforms to
[stats-contract](https://johanlassen.github.io/carpeDIAMS/reference/stats-contract.md);
pair with
[`apex_correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_filtering.md).

## See also

[`correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_statistics.md),
[`apex_correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_filtering.md)
