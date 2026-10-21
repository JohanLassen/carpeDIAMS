# Fast 1D intensity-adaptive grouping (IAAC), O(n log n)

Fast 1D intensity-adaptive grouping (IAAC), O(n log n)

## Usage

``` r
fast_1D_iaac_grouping_cpp(mz, intensity, k, p0 = 5, Isat = 50000)
```

## Arguments

- mz:

  numeric vector of m/z values (finite, \> 0)

- intensity:

  numeric vector (same length as mz; values \< 1 are clamped to 1)

- k:

  z-threshold; a group is kept iff every internal merge is \<= k sigma
  (k \< 0 yields all singletons; +Inf yields one group; NA is an error)

- p0:

  floor mass accuracy in ppm (high-intensity limit / systematic floor)

- Isat:

  intensity at which the per-point precision plateaus at p0

## Value

integer vector of group ids (length(mz)), 1-indexed by min original
position
