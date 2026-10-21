# Write a list of spectra as one MGF and optional SIRIUS .ms files

Standalone writer for spectra you assembled yourself; the pipeline uses
[`export_mgf()`](https://johanlassen.github.io/carpeDIAMS/reference/export_mgf.md).
Writes in chunks of `chunk_size` blocks, and one `feature_<id>.ms` file
per spectrum when `sirius_dir` is set.

## Usage

``` r
export_pseudospectra(
  spectra_list,
  mgf_path,
  sirius_dir = NULL,
  chunk_size = 1000L
)
```

## Arguments

- spectra_list:

  As described above.

- mgf_path:

  Output `.mgf` path.

- sirius_dir:

  Optional directory for per-feature `.ms` files. `NULL` (default) skips
  SIRIUS output entirely.

- chunk_size:

  Integer. MGF entries buffered before each `writeLines`.

## Value

Invisibly `NULL`.

## Details

Each element of `spectra_list` must be a list with:

- `peaks` — data.frame with columns `mz`, `total_intensity`

- `feature_id` — character, used as MGF `TITLE=` and SIRIUS `>compound`

- `mz_target` — numeric, precursor m/z

- `rtmed` — numeric, retention time in seconds

- `ESI` — character, "POSITIVE" / "NEGATIVE"
