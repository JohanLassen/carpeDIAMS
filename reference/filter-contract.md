# Filter contract

All filter closures in this package share one signature so they are
drop-in replaceable in
[`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md):

## Details

      function(stats) -> row-subset of stats

`stats` is the per-group table from a
[stats-contract](https://johanlassen.github.io/carpeDIAMS/reference/stats-contract.md)
closure. The filter returns the rows to keep, with the same columns.
Zero rows (or `NULL`) means the feature gets a
`DeconvERROR:NoAssociatedFragments` block in the MGF.
