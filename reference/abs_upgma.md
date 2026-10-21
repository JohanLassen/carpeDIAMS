# Absolute-tolerance UPGMA grouper

As
[`ppm_upgma()`](https://johanlassen.github.io/carpeDIAMS/reference/ppm_upgma.md),
but the maximum group width is a fixed value in Da. Uses only `mz`.

## Usage

``` r
abs_upgma(max_diff = 0.005)
```

## Arguments

- max_diff:

  Maximum group width (max m/z − min m/z) in Da.

## Value

A grouper closure - see
[grouper-contract](https://johanlassen.github.io/carpeDIAMS/reference/grouper-contract.md).

## Examples

``` r
g <- abs_upgma(max_diff = 0.005)
g(c(100.0, 100.001, 200.0))
#> [1] 1 1 3
```
