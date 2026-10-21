# Correlation + FDR + apex-agreement filter

As
[`correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_filtering.md),
and also requires `abs(apex_shift_med) <= max_apex_shift_med` and
`apex_shift_iqr <= max_apex_shift_iqr`. Errors if those columns are
missing, i.e. if
[`apex_correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_statistics.md)
was not used. Conforms to
[filter-contract](https://johanlassen.github.io/carpeDIAMS/reference/filter-contract.md).

## Usage

``` r
apex_correlation_filtering(
  cor_thresh = 0.3,
  fdr_thresh = 0.01,
  max_apex_shift_med = 0.5,
  max_apex_shift_iqr = 0.75
)
```

## Arguments

- cor_thresh:

  Minimum Pearson correlation to keep.

- fdr_thresh:

  Maximum BH-adjusted FDR to keep.

- max_apex_shift_med:

  Maximum `abs(apex_shift_med)`, in scans.

- max_apex_shift_iqr:

  Maximum `apex_shift_iqr`, in scans.

## Value

A filter closure — see
[filter-contract](https://johanlassen.github.io/carpeDIAMS/reference/filter-contract.md).

## See also

[`correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_filtering.md),
[`apex_correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_statistics.md)
