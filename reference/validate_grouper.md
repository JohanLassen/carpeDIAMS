# Sanity-check a custom grouper against the contract

Calls `grouper` on a small synthetic input and errors unless it returns
one id per input value and no more groups than distinct m/z values. It
does not check the type of the ids.

## Usage

``` r
validate_grouper(
  grouper,
  mz_example = c(100, 100.001, 100.5, 200),
  intensity = c(1000, 800, 200, 5000),
  target_mz = 100,
  scan = c(1L, 1L, 2L, 2L),
  ms1_int = c(5000, 5000, 4000, 4000)
)
```

## Arguments

- grouper:

  A function conforming to
  [grouper-contract](https://johanlassen.github.io/carpeDIAMS/reference/grouper-contract.md).

- mz_example:

  Numeric vector to probe with.

- intensity, target_mz, scan, ms1_int:

  Probe values passed to `grouper`.

## Value

Invisibly `TRUE`; errors otherwise.
