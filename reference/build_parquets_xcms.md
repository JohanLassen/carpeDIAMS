# Build pool + pointer parquets for one sample, from a subsetted msobject

The work behind `process_xcms(stage = "sample")`: looks up the sample's
peaks for `qualifying_features`, widens their RT windows by `rt_buffer`,
maps them to scans with
[`extractScanInfo()`](https://johanlassen.github.io/carpeDIAMS/reference/extractScanInfo.md),
and calls
[`Ms2PseudoRaw()`](https://johanlassen.github.io/carpeDIAMS/reference/Ms2PseudoRaw.md).
Writes nothing if the sample has no qualifying features or no matching
scans.

## Usage

``` r
build_parquets_xcms(
  mssubset,
  fd,
  out_dir,
  sample_id,
  qualifying_features,
  ppm = NULL,
  rt_buffer = 0
)
```

## Arguments

- mssubset:

  Output of
  [`subset_samples_xcms()`](https://johanlassen.github.io/carpeDIAMS/reference/subset_samples_xcms.md).

- fd:

  `featureDefinitions()` table.

- out_dir:

  Output directory (created if missing).

- sample_id:

  Index into `fileNames(mssubset$msobject)`.

- qualifying_features:

  Feature ids to build for this sample.

- ppm:

  MS1 tolerance in ppm, converted to Da at the largest feature m/z in
  this sample. `NULL` uses 0.005 Da.

- rt_buffer:

  Seconds added to both sides of each peak's RT window.

## Value

`NULL`.
