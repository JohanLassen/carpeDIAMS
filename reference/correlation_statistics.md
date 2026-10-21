# Per-group correlation statistics (the package default)

For each fragment group:

## Usage

``` r
correlation_statistics(min_count = 10L)
```

## Arguments

- min_count:

  Groups need more than this many `(sample, ms1_scan)` cells to be
  scored.

## Value

A statistics closure — see
[stats-contract](https://johanlassen.github.io/carpeDIAMS/reference/stats-contract.md).

## Details

1.  Sum fragment intensity per `(sample, ms1_scan)` cell.

2.  Drop the group unless it has **more than** `min_count` cells.

3.  Pearson correlation between the cell sums and the MS1 intensity,
    across all cells of all samples, with a two-sided t-test p-value.

4.  Benjamini–Hochberg FDR across the groups of this feature.

Output columns: `int` is the mean cell sum, `count` the number of cells,
`n_samples` the number of samples. Conforms to
[stats-contract](https://johanlassen.github.io/carpeDIAMS/reference/stats-contract.md);
pair with
[`correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_filtering.md).

## See also

[`apex_correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_statistics.md),
[`correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_filtering.md)
