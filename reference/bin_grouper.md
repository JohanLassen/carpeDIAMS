# Fixed-width binning grouper

Group id is `floor(mz / bin_width) + 1`. Fast, but two close values on
either side of a bin edge end up in different groups. Uses only `mz`.

## Usage

``` r
bin_grouper(bin_width = 0.01)
```

## Arguments

- bin_width:

  Bin width in Da.

## Value

A grouper closure - see
[grouper-contract](https://johanlassen.github.io/carpeDIAMS/reference/grouper-contract.md).

## Examples

``` r
g <- bin_grouper(bin_width = 0.01)
g(c(100.001, 100.009, 100.015))
#> [1] 10001 10001 10002
```
