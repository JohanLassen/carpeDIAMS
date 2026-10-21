# ppm-based UPGMA grouper (the package default)

1-D UPGMA on the sorted m/z values, cut so that no group is wider (max
m/z − min m/z) than `ppm` of the precursor m/z. Uses `mz` and
`target_mz`.

## Usage

``` r
ppm_upgma(ppm = 30)
```

## Arguments

- ppm:

  Maximum group width, in ppm of `target_mz`.

## Value

A grouper closure - see
[grouper-contract](https://johanlassen.github.io/carpeDIAMS/reference/grouper-contract.md).

## See also

[`abs_upgma()`](https://johanlassen.github.io/carpeDIAMS/reference/abs_upgma.md),
[`iaac_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/iaac_grouper.md),
[`pelt_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/pelt_grouper.md),
[`validate_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/validate_grouper.md)

## Examples

``` r
g <- ppm_upgma(ppm = 30)
g(c(100.0, 100.001, 200.0), target_mz = 100)
#> [1] 1 1 3
```
