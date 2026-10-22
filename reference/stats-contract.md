# Statistics contract

All group-statistics closures in this package share one signature so
they are drop-in replaceable in
[`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md)
and
[`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md):

## Details

      function(dt) -> per-group summary data.table

`dt` holds one feature's rows from
[`batch_get_feature_data()`](https://johanlassen.github.io/carpeDIAMS/reference/batch_get_feature_data.md)
(one row per MS1 scan × following MS2 centroid), after
[`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md)
has dropped rows without a fragment, applied the apex window, and added
`group` (from the grouper) and `group_mz` (the consensus m/z, constant
within a group).

The closure returns one row per group with at least

      group, fragment_mz, mz, int, ms1_int,
      correlation, p_value, fdr, count, n_samples

where `fragment_mz` is `group_mz` carried through. Extra columns are
allowed (e.g. `apex_shift_med`, `apex_shift_iqr`). It should not
threshold; that is the
[filter-contract](https://johanlassen.github.io/carpeDIAMS/reference/filter-contract.md)
closure's job.
[`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md)
reads `fragment_mz`, `group`, `int`, `n_samples`, `count`,
`correlation`, `p_value` and `fdr` from the filtered result.

Optionally attach the `(group, sample, ms1_scan)` aggregate as attribute
`"grouped_signal"`; the diagnostic report plots it.
