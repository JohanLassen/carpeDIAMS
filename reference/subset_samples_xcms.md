# Pick representative samples with apricot MaxCoverage

A sample "covers" a feature when its intensity is at or above that
feature's `q_prob` quantile across samples (missing values count as 0).
apricot's `MaxCoverageSelection` then picks the `subset_size` samples
that together cover the most features. Needs Python with
`apricot-select` and `numpy` (via reticulate).

## Usage

``` r
subset_samples_xcms(
  msobject,
  msassay,
  subset_size = 100L,
  q_prob = 0.9,
  subsample = TRUE
)
```

## Arguments

- msobject:

  An `XcmsExperiment`/`XCMSnExp`.

- msassay:

  Samples × features intensity matrix; `rownames(msassay)` must be
  `basename(xcms::fileNames(msobject))`.

- subset_size:

  Number of samples to keep.

- q_prob:

  Per-feature quantile that defines "covered".

- subsample:

  If `FALSE`, return every sample and skip Python.

## Value

List with `msobject` (subsetted), `assay` (selected rows), `peaks` and
`all_files` (both from the **unsubsetted** object, so `peakidx` and
sample indices stay valid).
