# Exact 1D grouping via PELT (penalized optimal partitioning)

Exact 1D grouping via PELT (penalized optimal partitioning)

## Usage

``` r
fast_1D_pelt_grouping_cpp(
  mz,
  intensity,
  beta = -1,
  p0 = 5,
  Isat = 50000,
  max_seg = 0L
)
```

## Arguments

- mz:

  numeric vector of m/z values (finite, \> 0)

- intensity:

  numeric vector (same length as mz; values \< 1 are clamped to 1)

- beta:

  penalty per additional segment; beta \<= 0 -\> auto (2 \* log(n))

- p0:

  floor mass accuracy in ppm (high-intensity limit / systematic floor)

- Isat:

  intensity at which the per-point precision plateaus at p0

- max_seg:

  maximum number of centroids allowed in one group; \<= 0 = unlimited.
  Caps PELT's worst case (one giant group) at O(n \* max_seg) instead of
  O(n^2). The result stays the exact global optimum over partitions
  obeying the cap.

## Value

integer vector of group ids (length(mz)), 1-indexed by min original
position
