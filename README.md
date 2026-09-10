
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- knitr::knit("README.Rmd") -->

# carpeDIAMS

<!-- badges: start -->
<!-- badges: end -->

**carpeDIAMS** turns DIA (bbCID / dia-PASEF-style) mass spectrometry
acquisitions into pseudo-DDA spectra. It correlates MS1 precursor
elution with MS2 fragment elution across samples to figure out which
fragments actually belong to which precursor — something DIA doesn’t
tell you directly, because every MS2 scan fragments *everything*
co-isolated in the window, not just your compound of interest.

The output is a standard MGF: one block per feature, ready for SIRIUS,
GNPS, or any downstream MS2-based annotation tool.

## Installation

``` r
pak::pak("JohanLassen/carpeDIAMS")
```

`carpeDIAMS` also needs `apricot-select` + `numpy` on the Python side if
you plan to use `subsample = TRUE` (large-cohort sample selection, see
below) — everything else is pure R/Rcpp and works without Python.

## The three inputs, one function

Whatever stage your data is at, one function call gets you to an MGF.
`carpediams()` looks at what you hand it and routes to the matching
worker:

| You have                                            | Route                  |
|:----------------------------------------------------|:-----------------------|
| A single mzML, exploratory, no prior knowledge      | `process_mzml()`       |
| mzML(s) + a target list of known compounds (mz, rt) | `process_targetlist()` |
| A grouped `XcmsExperiment` / `XCMSnExp`             | `process_xcms()`       |

You can call the specific function directly, or just call `carpediams()`
and let it figure out which one you meant:

``` r
library(carpeDIAMS)

# one exploratory file, no target list -> process_mzml()
carpediams(
  input  = "sample.mzML",
  output = "sample.mgf"
)

# mzML(s) against a target list of known mz/rt -> process_targetlist()
carpediams(
  input      = c("sample_a.mzML", "sample_b.mzML"),
  output     = "targets.mgf",
  targetlist = "targets.csv"        # needs `mz` and `rt` columns
)

# an already-grouped xcms object -> process_xcms()
carpediams(
  input  = "grouped_features.rds",  # or an in-memory XcmsExperiment
  output = "features.mgf"
)
```

That’s the whole interface for everyday use. Everything below is about
what happens under the hood and how to change it.

## What each route does

**`process_mzml()`** — loads one mzML, detects chromatographic peaks
with centWave, and builds pseudo-spectra for every peak found. No prior
knowledge required; slowest to reach an answer since it does peak
detection from scratch, but needs nothing but the raw file.

**`process_targetlist()`** — skips peak detection entirely. Given a CSV
of `mz`/`rt` pairs (with an optional `feature` id column), it extracts
scans in a window around each target’s retention time and builds
pseudo-spectra only for those. Fast and lightweight — ideal when you
already know what you’re looking for (a drug panel, a reference list of
standards).

**`process_xcms()`** — the multi-sample workhorse. Takes an
`XcmsExperiment`/`XCMSnExp` that has already been through
`xcms::groupChromPeaks()`, and builds pseudo-spectra for every feature,
pooling evidence across every sample that feature appears in. This is
what you want once you’re past exploration and running a real
multi-sample study.

All three converge on the same underlying machinery — extraction,
grouping, statistics, filtering, MGF export — so the tuning knobs below
work identically regardless of which route you took.

## Why pseudo-spectra need statistics at all

In bbCID/DIA, an MS2 scan doesn’t isolate one precursor — it fragments
everything eluting at that moment. So the pool of MS2 fragments near
your feature’s retention time is a mixture: real fragments of your
compound, plus fragments from anything else that happened to coelute.
carpeDIAMS’s job is separating those two populations, per feature, per
sample, using the shape of the MS1 elution as ground truth: a true
fragment’s intensity across scans should track the precursor’s MS1
intensity across scans.

Concretely, for each feature:

1.  **Extract** raw MS1 + MS2 rows around the feature’s expected
    retention time.
2.  **Group** MS2 centroids into fragment candidates (`grouper`) — the
    the m/z clustering step.
3.  **Score** each candidate group with per-group statistics
    (`group_statistics_fn`) — typically the Pearson correlation between
    summed MS2 intensity and MS1 intensity across scans, with a
    Benjamini-Hochberg FDR correction.
4.  **Filter** to the surviving groups (`filter_fn`) — threshold on
    those statistics.
5.  **Export** surviving groups as one MGF block per feature.

Steps 2-4 are each a swappable closure. That’s deliberate: which
grouping/statistics/filtering approach is *right* is an open research
question, not a solved one, so the package treats it as a pluggable
layer rather than baking in one fixed algorithm.

## Tuning the pipeline

### Grouping (`grouper`)

Controls how raw MS2 centroids get clustered into fragment candidates.

``` r
ppm_upgma(ppm = 30)                       # default: ppm-scaled UPGMA
abs_upgma(max_diff = 0.005)               # fixed-Da UPGMA
iaac_grouper(k = 1, p0 = 8, Isat = 1e5)   # intensity-adaptive clustering
pelt_grouper(beta = -1, p0 = 5, Isat = 5e4)  # globally-optimal PELT partition
bin_grouper(bin_width = 0.01)             # fixed-width binning (fast, crude)
```

``` r
carpediams(input = "grouped.rds", output = "out.mgf",
           grouper = iaac_grouper(k = 1, p0 = 8, Isat = 1e5))
```

### Statistics + filtering (`group_statistics_fn`, `filter_fn`)

Two closures, always used as a pair. The statistics closure computes
per-group numbers; the filter closure thresholds on them. They’re
separate so you can compute richer statistics without committing to one
fixed cutoff logic.

``` r
correlation_statistics(min_count = 10L)      # default: Pearson r + BH FDR
apex_correlation_statistics(min_count = 10L) # + intra-sample MS1<->MS2 apex agreement

correlation_filtering(cor_thresh = 0.3, fdr_thresh = 0.01)         # default
apex_correlation_filtering(cor_thresh = 0.3, fdr_thresh = 0.01,
                           max_apex_shift_med = 0.5,
                           max_apex_shift_iqr = 0.75)
```

The apex-agreement pair is worth knowing about: it measures, per sample,
whether the intensity-weighted MS1 apex and MS2 apex line up in
scan-index space. A real fragment tracks its precursor’s elution apex; a
coeluting contaminant usually doesn’t, even when it happens to correlate
well overall. Use it alongside (or instead of) plain correlation when
correlation alone lets through fragments from coeluting compounds:

``` r
carpediams(
  input  = "grouped.rds",
  output = "out.mgf",
  group_statistics_fn = apex_correlation_statistics(),
  filter_fn           = apex_correlation_filtering(max_apex_shift_med = 0.5,
                                                    max_apex_shift_iqr = 0.75)
)
```

`apex_correlation_filtering()` requires the two apex columns to be
present — pair it with `apex_correlation_statistics()`, not the plain
`correlation_statistics()`. If you forget, you get a clear error telling
you which columns are missing and what to add.

### Writing your own

All three closures are plain functions with a documented contract — see
`?grouper-contract`, `?stats-contract`, and `?filter-contract`.
Implementing your own is the expected way to try a new idea:

``` r
my_filter <- function(stats) {
  stats[correlation > 0.6 & count >= 20]
}

carpediams(input = "grouped.rds", output = "out.mgf",
           filter_fn = my_filter)
```

`validate_grouper()` sanity-checks a custom grouper against the contract
before you wire it into a full run.

## Scaling up

Everything above runs single-process, single-machine. When a cohort gets
large enough that a single `xcms::groupChromPeaks()` run doesn’t
comfortably fit one job, `process_xcms()` (and `process_targetlist()`)
expose a `stage` argument that splits the same pipeline into
DAG-schedulable pieces — one job per sample, batched filtering across
many jobs, with a documented on-disk contract in between. See
`vignette("dag-cli", package = "carpeDIAMS")` for the full walkthrough,
including a Snakemake sketch and a failure/resume playbook.

If you have thousands of samples and don’t need every one of them
represented per feature,
`process_xcms(..., subsample = TRUE, subset_size = 100)` uses apricot’s
`MaxCoverageSelection` to pick a representative subset before building
pseudo-spectra, instead of paying the cost of every sample.

## Learn more

- `vignette("carpeDIAMS")` — this document, as an installed vignette.
- `vignette("dag-cli")` — running carpeDIAMS across a cluster.
- `?carpediams` — full argument reference and routing rules.
- `?clean_pseudo_spectra` — the batch orchestrator behind every route,
  if you need lower-level control (custom DuckDB connections, MGF
  sidecar output per feature, PDF diagnostic reports).
