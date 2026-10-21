# Legacy per-group thresholder (superseded by [`correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_filtering.md))

Not exported. Use `correlation_filtering(cor_thresh, fdr_thresh)`.

## Usage

``` r
select_pseudo_spec(correlations, cor_thresh = 0.5, fdr_thresh = 0.05)
```

## Arguments

- correlations:

  Per-group summary from a
  [stats-contract](https://johanlassen.github.io/carpeDIAMS/reference/stats-contract.md)
  closure (typically
  [`correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_statistics.md)).

- cor_thresh:

  Numeric. Minimum correlation to keep.

- fdr_thresh:

  Numeric. Maximum FDR to keep.

## Value

A `data.table` (possibly empty).
