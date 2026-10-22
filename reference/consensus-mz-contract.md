# Consensus m/z contract

All consensus-m/z closures in this package share one signature so they
are drop-in replaceable in
[`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md):

## Details

      function(mz, int) -> single numeric consensus m/z

Called once per fragment group, after grouping; the result becomes the
group's `fragment_mz`. Should return `NA_real_`, not error, when the
group has no finite values.
