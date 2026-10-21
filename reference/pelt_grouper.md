# PELT (penalized optimal partitioning) grouper

Exact optimal partition of the sorted m/z values: minimises the summed
per-group cost plus `beta` per group. Per-centroid uncertainty is the
same as in
[`iaac_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/iaac_grouper.md).
Larger `beta` gives fewer, wider groups. Uses `mz` and `intensity`. See
[`fast_1D_pelt_grouping_cpp()`](https://johanlassen.github.io/carpeDIAMS/reference/fast_1D_pelt_grouping_cpp.md).

## Usage

``` r
pelt_grouper(beta = -1, p0 = 5, Isat = 50000, max_seg = 0L)
```

## Arguments

- beta:

  Penalty per group; `<= 0` uses `2 * log(n)`.

- p0:

  Mass accuracy floor (ppm) reached at high intensity.

- Isat:

  Intensity at which the uncertainty reaches `p0`.

- max_seg:

  Maximum centroids per group (`0` = no limit). Bounds run time at
  `O(n * max_seg)`; the result is optimal under the cap.

## Value

A grouper closure - see
[grouper-contract](https://johanlassen.github.io/carpeDIAMS/reference/grouper-contract.md).

## See also

[`iaac_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/iaac_grouper.md),
[`ppm_upgma()`](https://johanlassen.github.io/carpeDIAMS/reference/ppm_upgma.md),
[`abs_upgma()`](https://johanlassen.github.io/carpeDIAMS/reference/abs_upgma.md),
[`validate_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/validate_grouper.md)

## Examples

``` r
g <- pelt_grouper(beta = -1, p0 = 5, Isat = 5e4)
g(c(100.000, 100.001, 100.500), intensity = c(1000, 800, 200))
#> [1] 1 1 3
```
