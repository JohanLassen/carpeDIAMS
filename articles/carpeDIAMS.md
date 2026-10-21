# Getting started with carpeDIAMS

``` r

library(carpeDIAMS)
```

## Inputs

`carpediams(input, output, targetlist = NULL, ...)` picks a route from
its input:

| `input` | `targetlist` | Route |
|:---|:---|:---|
| `XcmsExperiment`/`XCMSnExp`, or `.rds`/`.qs`/`.rda` path | ignored | [`process_xcms()`](https://johanlassen.github.io/carpeDIAMS/reference/process_xcms.md) |
| mzML file, mzML directory, or vector of paths | given | [`process_targetlist()`](https://johanlassen.github.io/carpeDIAMS/reference/process_targetlist.md) |
| one mzML file | `NULL` | [`process_mzml()`](https://johanlassen.github.io/carpeDIAMS/reference/process_mzml.md) |
| mzML directory or vector of paths | `NULL` | error |

Any other argument goes to the routed function, so you can call that
function directly instead if you prefer.

### xcms object

``` r

carpediams(input = "grouped_xcms.rds", output = "features.mgf",
           top_n = 10, ppm = 20, rt_buffer = 5)
```

The object must have been through
[`xcms::groupChromPeaks()`](https://rdrr.io/pkg/xcms/man/groupChromPeaks.html).
Each feature in `featureDefinitions()` is built from its `top_n` most
intense samples (by `featureValues(value = "into")`). In each of those
samples, the scans come from the feature’s chromatographic peak, widened
by `rt_buffer` seconds. For large cohorts,
`subsample = TRUE, subset_size = 100` first reduces the samples to a
representative set with apricot (Python).

### mzML files and a target list

``` r

carpediams(input = "data/mzml/", output = "targets.mgf",
           targetlist = "targets.csv", rt_width = 15, ppm = 20)
```

The target list needs `mz` and `rt` columns, and optionally `feature`
for ids. An `rt` column whose maximum is below 60 is read as minutes.
Every target is extracted from every file, in an `rt_width`-second
window centred on its `rt`. No peak detection is done.

### One mzML file

``` r

carpediams(input = "sample.mzML", output = "sample.mgf")
```

[`process_mzml()`](https://johanlassen.github.io/carpeDIAMS/reference/process_mzml.md)
keeps 30–1200 s (`rt_range`), runs centWave with fixed parameters, and
builds a spectrum for every peak found. Fragments are scored across the
scans of that single file. This route accepts only `rt_range`,
`grouper`, `group_statistics_fn` and `filter_fn`.

## Defaults

| Step | Default | Change with |
|:---|:---|:---|
| Scan window per sample | apex = highest MS1 scan; keep scans \< 50 scans from it | `buffer`, `apex_source`, `window_mode` |
| m/z grouping | `ppm_upgma(ppm = 20)` (xcms, target list) or `ppm_upgma(ppm = 30)` (mzML) | `grouper`, `ppm` |
| Consensus m/z | `consensus_mz_weighted(frac = 0.5, Isat = 1e4)` | `consensus_mz_fn` |
| Statistics | `correlation_statistics(min_count = 10)` | `group_statistics_fn`, `min_count` |
| Filter | `correlation_filtering(cor_thresh = 0.1, fdr_thresh = 0.1)` | `filter_fn`, `cor_thresh`, `fdr_thresh` |

`buffer`, `apex_source`, `window_mode` and `consensus_mz_fn` are passed
through `...` and work on the xcms and target-list routes only.

The pipeline’s default filter (0.1 / 0.1) is more permissive than the
defaults of the
[`correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_filtering.md)
factory itself (0.3 / 0.01). If you build the filter yourself, set both
thresholds explicitly.

## Choosing a grouper

All samples’ centroids for a feature are grouped in one call, so the
tolerance has to cover both scatter within a sample and calibration
drift between samples.

``` r

ppm_upgma(ppm = 20)                       # group width <= 20 ppm of the precursor m/z
abs_upgma(max_diff = 0.005)               # group width <= 0.005 Da
iaac_grouper(k = 2, p0 = 12, Isat = 1e4)  # intensity-dependent tolerance
pelt_grouper(p0 = 15, Isat = 1e5)         # optimal partition, same noise model
bin_grouper(bin_width = 0.01)             # fixed bins
```

[`iaac_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/iaac_grouper.md)
and
[`pelt_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/pelt_grouper.md)
give each centroid an m/z uncertainty of `p0` ppm at high intensity,
growing as `sqrt(Isat / intensity)` below `Isat`. Weak centroids are
therefore allowed to scatter more.

## Statistics and filters

For each fragment group,
[`correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_statistics.md)
sums the fragment intensity per (sample, MS1 scan), correlates those
sums with the MS1 intensity, and reports `correlation`, `p_value`,
BH-adjusted `fdr`, `count` (number of cells) and `n_samples`. Groups
need more than `min_count` cells.

[`apex_correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_statistics.md)
also reports how far the fragment’s elution centre is from the
precursor’s, in scans (`apex_shift_med`, `apex_shift_iqr` across
samples). Pair it with
[`apex_correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_filtering.md)
to reject co-eluting fragments that correlate but peak at the wrong
time:

``` r

carpediams(
  input = "grouped_xcms.rds", output = "features.mgf",
  group_statistics_fn = apex_correlation_statistics(min_count = 10),
  filter_fn = apex_correlation_filtering(cor_thresh = 0.5, fdr_thresh = 0.01,
                                         max_apex_shift_med = 0.5,
                                         max_apex_shift_iqr = 0.75)
)
```

A filter is just a function of the statistics table:

``` r

carpediams(input = "grouped_xcms.rds", output = "features.mgf",
           filter_fn = function(stats) stats[correlation > 0.6 & count >= 20])
```

## The MGF

Every feature (or target) gets exactly one block, in input order, with
`TITLE`, `PEPMASS`, `RTINSECONDS`, `CHARGE=1+`, `MSLEVEL=2` and
`IONMODE`. Peak intensity is the mean, over (sample, scan) cells, of the
summed fragment intensity. When no spectrum could be built, the block
has no peaks and the title says why:

| Title suffix | Meaning |
|:---|:---|
| `DeconvERROR:MissingTargetInfo` | no precursor m/z (or no RT on the target-list route) |
| `DeconvERROR:MissingMs1Signal` | no MS1 scans were extracted for the feature |
| `DeconvERROR:NoAssociatedFragments` | no fragment group passed the statistics and filter |

## Larger datasets

[`carpediams()`](https://johanlassen.github.io/carpeDIAMS/reference/carpediams.md)
runs everything in one R session by default. The xcms and target-list
routes can also run as separate stages (`prepare`, one `sample` job per
file, `filter` in batches, `finalize`) that share a working directory:

- [`vignette("xcms-interactive")`](https://johanlassen.github.io/carpeDIAMS/articles/xcms-interactive.md)
  walks through the stages in R.
- [`vignette("dag-cli")`](https://johanlassen.github.io/carpeDIAMS/articles/dag-cli.md)
  covers running them as cluster jobs.
