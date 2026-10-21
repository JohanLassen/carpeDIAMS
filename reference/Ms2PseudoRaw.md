# Build pool + pointer parquets for one sample

Reads the scans of one file between the lowest and highest scan in
`single_sample` (± `buffer`) and writes two parquets to `output_dir`:

## Usage

``` r
Ms2PseudoRaw(
  single_sample,
  spectra_object,
  output_dir = "results/",
  maxDiff = 0.01,
  buffer = 50L,
  return_spectra = FALSE
)
```

## Arguments

- single_sample:

  One file's rows from
  [`extractScanInfo()`](https://johanlassen.github.io/carpeDIAMS/reference/extractScanInfo.md)
  (or the same columns built from a target list).

- spectra_object:

  `Spectra` object containing that file.

- output_dir:

  Output directory.

- maxDiff:

  MS1 m/z tolerance in Da.

- buffer:

  Extra scans on each side.

- return_spectra:

  If `TRUE`, also return the two tables.

## Value

`NULL`, or `list(pool, pointer)` when `return_spectra = TRUE`. Returns
`NULL` without writing if no scans are in range.

## Details

- `sample_<id>_pool.parquet`: every MS2 centroid in that scan range
  (`sample, scan, rt, mz, int`).

- `sample_<id>_pointers.parquet`: for each feature and each MS1 scan
  within its scan range ± `buffer`, the summed intensity of MS1
  centroids within `maxDiff` Da of the precursor
  (`feature, scan, mz_target, apex_scan, ms1_int, mz, rt, sample`).
