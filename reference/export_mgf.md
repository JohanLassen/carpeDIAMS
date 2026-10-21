# Assemble a combined MGF from a flat pseudo-spectra results table

Writes one MGF block per row of `fd`, in `fd` order. Features with
fragment rows get `TITLE`, `PEPMASS` (from `mzmed`), `RTINSECONDS` (from
`rtmed`), `CHARGE=1+`, `MSLEVEL=2`, `IONMODE` and the peaks sorted by
m/z. Features with an error row, or missing from `results`, get an empty
block titled `<feature_id> DeconvERROR:<reason>`
(`NoAssociatedFragments` when missing). Warns if no feature has peaks.

## Usage

``` r
export_mgf(results, fd, mgf_path, ion_mode = "POSITIVE")
```

## Arguments

- results:

  Output of
  [`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md):
  `feature_id`, `mz`, `total_intensity`, optional `deconv_error`.

- fd:

  Feature table with `mzmed`/`rtmed` (or `mz`/`rt`); rownames or a
  `feature_id` column give the ids.

- mgf_path:

  Output path (directories are created).

- ion_mode:

  Written as `IONMODE=`.

## Value

Invisibly `mgf_path`.
