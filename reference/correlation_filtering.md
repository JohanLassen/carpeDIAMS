# Correlation + BH FDR filter (the package default)

Keeps groups with `correlation >= cor_thresh` and `fdr <= fdr_thresh`.
Groups with `NA` correlation are dropped. Conforms to
[filter-contract](https://johanlassen.github.io/carpeDIAMS/reference/filter-contract.md).

## Usage

``` r
correlation_filtering(cor_thresh = 0.3, fdr_thresh = 0.01)
```

## Arguments

- cor_thresh:

  Minimum Pearson correlation to keep.

- fdr_thresh:

  Maximum BH-adjusted FDR to keep.

## Value

A filter closure — see
[filter-contract](https://johanlassen.github.io/carpeDIAMS/reference/filter-contract.md).

## Details

The factory defaults (0.3, 0.01) differ from the defaults the
`process_*()` functions build when `filter_fn = NULL` (0.1, 0.1).

## See also

[`apex_correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_filtering.md),
[`correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_statistics.md)
