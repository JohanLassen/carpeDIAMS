# Intensity-adaptive agglomerative clustering grouper

Agglomerative 1-D clustering where each centroid's m/z uncertainty
depends on its intensity:
`sigma_i = mz_i * p0 * max(1, sqrt(Isat / I_i)) / 1e6`. A group is kept
only if every merge inside it is within `k` sigma. Uses `mz` and
`intensity`. See
[`fast_1D_iaac_grouping_cpp()`](https://johanlassen.github.io/carpeDIAMS/reference/fast_1D_iaac_grouping_cpp.md).

## Usage

``` r
iaac_grouper(k = 1, p0 = 8, Isat = 1e+05)
```

## Arguments

- k:

  Merge threshold in sigma. `k < 0` gives all singletons, `Inf` one
  group.

- p0:

  Mass accuracy floor (ppm) reached at high intensity.

- Isat:

  Intensity at which the uncertainty reaches `p0`.

## Value

A grouper closure - see
[grouper-contract](https://johanlassen.github.io/carpeDIAMS/reference/grouper-contract.md).

## See also

[`ppm_upgma()`](https://johanlassen.github.io/carpeDIAMS/reference/ppm_upgma.md),
[`abs_upgma()`](https://johanlassen.github.io/carpeDIAMS/reference/abs_upgma.md),
[`pelt_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/pelt_grouper.md),
[`validate_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/validate_grouper.md)

## Examples

``` r
g <- iaac_grouper(k = 1, p0 = 8, Isat = 1e5)
g(c(100.0, 100.001, 200.0), intensity = c(1000, 800, 200))
#> [1] 1 1 3
```
