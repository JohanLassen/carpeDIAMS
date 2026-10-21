# Write a single-feature MGF block (concat-friendly)

Writes one `BEGIN IONS ... END IONS` block (`TITLE`, `PEPMASS`, optional
`RTINSECONDS`, `CHARGE=1+`, peaks). Used by
[`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md)
when `mgf_dir` is set. Files can be joined with `cat`.

## Usage

``` r
write_simple_mgf(pseudo_spec, mz_target, feature_id, path, rt = NA_real_)
```

## Arguments

- pseudo_spec:

  Filtered per-group table with `fragment_mz` and `int`.

- mz_target:

  Numeric. Precursor m/z.

- feature_id:

  Used as `TITLE=`.

- path:

  Output `.mgf` path.

- rt:

  Optional retention time (seconds); skipped when `NA`.

## Value

Invisibly `path`, or `NULL` when `pseudo_spec` is empty.
