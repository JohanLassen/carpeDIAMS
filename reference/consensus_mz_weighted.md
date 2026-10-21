# Intensity-weighted FWHM-core consensus m/z (the package default)

Intensity-weighted mean m/z of the centroids with intensity
`>= frac * max(intensity)` in the group, with weights capped at `Isat`.
Low-intensity centroids are dropped because their m/z is least precise.
Falls back to the unweighted mean of those centroids if all weights are
zero.

## Usage

``` r
consensus_mz_weighted(frac = 0.5, Isat = 10000)
```

## Arguments

- frac:

  Fraction of the group's maximum intensity a centroid needs to be
  included (`0.5` = half maximum).

- Isat:

  Weight cap, so a few saturated centroids cannot dominate.

## Value

A consensus-mz closure - see
[consensus-mz-contract](https://johanlassen.github.io/carpeDIAMS/reference/consensus-mz-contract.md).

## See also

[`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md)

## Examples

``` r
f <- consensus_mz_weighted(frac = 0.5, Isat = 1e4)
f(c(100.000, 100.001, 100.500), c(1000, 800, 50))
#> [1] 100.0004
```
