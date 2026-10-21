# Grouper contract

All groupers in this package share one signature so they are drop-in
replaceable in
[`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md)
and
[`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md):

## Details


      function(mz, intensity = NULL, target_mz = NULL, scan = NULL, ms1_int = NULL)

Only `mz` is required; each grouper ignores the arguments it does not
use.
[`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md)
calls the grouper **once per feature**, on the MS2 centroids of all
samples pooled together, passing `mz` (fragment m/z), `intensity`
(centroid intensity) and `target_mz` (precursor m/z). `scan` and
`ms1_int` are not passed there; they exist for custom callers.

Because samples are pooled, the tolerance has to cover between-sample
calibration drift as well as within-sample scatter.

A grouper returns an integer vector, one group id per element of `mz`,
with no more distinct ids than distinct `mz` values. Check a custom
grouper with
[`validate_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/validate_grouper.md).
