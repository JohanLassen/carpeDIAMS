# carpeDIAMS

carpeDIAMS builds MS2 spectra for LC-MS features acquired with broadband
CID (bbCID), where MS1 and all-ion MS2 scans alternate and no precursor
is isolated. Every MS2 scan contains fragments of everything that
co-elutes, so carpeDIAMS keeps a fragment only if its intensity rises
and falls with the precursor’s MS1 signal. The result is an MGF file
with one pseudo-DDA spectrum per feature, ready for spectral library
search, SIRIUS or molecular networking.

## Installation

``` r
# install.packages("pak")
pak::pak("JohanLassen/carpeDIAMS")
```

Sample subsampling (`subsample = TRUE`) also needs Python with
`apricot-select` and `numpy`, found through reticulate.

## Usage

[`carpediams()`](https://johanlassen.github.io/carpeDIAMS/reference/carpediams.md)
takes one of three kinds of input and writes an MGF:

``` r
library(carpeDIAMS)

# A grouped xcms object (XcmsExperiment or .rds/.qs/.rda path)
carpediams(input = "grouped_xcms.rds", output = "features.mgf")

# mzML files plus a target list with mz and rt columns
carpediams(input = c("a.mzML", "b.mzML"), output = "targets.mgf",
           targetlist = "targets.csv")

# A single mzML file (centWave peak detection is run for you)
carpediams(input = "sample.mzML", output = "sample.mgf")
```

## How a spectrum is built

For each feature:

1.  For each sample, take the MS1 scans around the feature’s peak and
    the MS2 scan that follows each one. Record the precursor’s MS1
    intensity and the MS2 centroids up to precursor m/z + 2.
2.  Group the MS2 centroids of all samples by m/z (`grouper`).
3.  Per group, correlate the summed fragment intensity with the MS1
    intensity across scans and samples (`group_statistics_fn`).
4.  Keep groups that pass a correlation and FDR threshold (`filter_fn`).

Steps 2–4 are plain R functions you can replace. Features where nothing
survives still get an MGF entry, titled with a `DeconvERROR` reason.

## Learn more

- [`vignette("carpeDIAMS")`](https://johanlassen.github.io/carpeDIAMS/articles/carpeDIAMS.md):
  inputs, defaults and tuning.
- [`vignette("xcms-interactive")`](https://johanlassen.github.io/carpeDIAMS/articles/xcms-interactive.md):
  the staged xcms workflow, step by step.
- [`vignette("dag-cli")`](https://johanlassen.github.io/carpeDIAMS/articles/dag-cli.md):
  running the stages as cluster jobs.
- [`vignette("developer")`](https://johanlassen.github.io/carpeDIAMS/articles/developer.md):
  comparing and writing groupers and filters.
