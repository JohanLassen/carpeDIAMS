# Evaluating a custom grouper (developer guide)

Which grouping/filtering approach is *right* for a given instrument,
acquisition mode, or compound class is an open research question, not
something carpeDIAMS tries to settle once and for all. The `grouper`,
`group_statistics_fn`, and `filter_fn` arguments to
\[[`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md)\]/\[[`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md)\]
exist so you can swap in your own answer without forking the package.
This vignette walks through evaluating one, using a real, shipped
example rather than a synthetic one.

## The fixture

`inst/testdata/FT00964/` is a real slice of a production carpeDIAMS run
— not synthetic data. It’s feature `FT00964` (precursor m/z ≈ 166.05),
across the 10 samples that carry it, trimmed to exactly the rows the
real query in
[`batch_get_feature_data()`](https://johanlassen.github.io/carpeDIAMS/reference/batch_get_feature_data.md)
ever touches for this feature (pool/pointers parquets, same on-disk
contract as any other carpeDIAMS run — see
[`vignette("dag-cli")`](https://johanlassen.github.io/carpeDIAMS/articles/dag-cli.md)
for what that contract looks like at scale). It’s small (~600 KB)
specifically so it can ship with the package and run in this vignette.

``` r

library(carpeDIAMS)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%

fixture_dir <- system.file("testdata", "FT00964", package = "carpeDIAMS")
fd <- readRDS(file.path(fixture_dir, "fd.rds"))
fd[, c("mzmed", "rtmed")]
#>            mzmed    rtmed
#> FT00964 166.0499 172.0318
```

Pull the raw MS1+MS2 rows the same way
[`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md)
does internally — connect, then fetch:

``` r

con <- setup_connection(fixture_dir)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpehZATr/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
raw <- batch_get_feature_data("FT00964", con)
DBI::dbDisconnect(con, shutdown = TRUE)

nrow(raw)
#> [1] 31455
length(unique(raw$sample))
#> [1] 10
```

`raw` has one row per (sample, MS1 scan, matched MS2 fragment) — exactly
what
\[[`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md)\]
expects as `dt`.

## Comparing groupers

Try the default (`ppm_upgma`) against `iaac_grouper` at the production
tolerance:

``` r

mz_target <- fd["FT00964", "mzmed"]

stats_ppm  <- process_pseudo_spec(raw, mz_target = mz_target,
                                  grouper = ppm_upgma(ppm = 30))
stats_iaac <- process_pseudo_spec(raw, mz_target = mz_target,
                                  grouper = iaac_grouper(k = 2, p0 = 12, Isat = 1e4))

# top fragment by correlation, under each grouper
stats_ppm[order(-abs(correlation))][1, .(fragment_mz, count, correlation)]
#>    fragment_mz count correlation
#>          <num> <int>       <num>
#> 1:    67.05413    96    0.916842
stats_iaac[order(-abs(correlation))][1, .(fragment_mz, count, correlation)]
#>    fragment_mz count correlation
#>          <num> <int>       <num>
#> 1:    67.05413    96    0.916842
```

Both agree closely here. That’s not always true, and this fixture is a
good place to see exactly why.

## A concrete `p0` sensitivity example

[`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md)
groups the centroids of all samples in one call, so the tolerance has to
absorb between-sample calibration drift as well as within-sample scan
noise. `p0` (for `iaac_grouper`/`pelt_grouper`) sets that tolerance. Too
tight, and a single real fragment gets reported as two:

``` r

for (p0 in c(4, 8, 12)) {
  g <- iaac_grouper(k = 1, p0 = p0, Isat = 1e5)
  s <- process_pseudo_spec(raw, mz_target = mz_target, grouper = g)
  near67 <- s[fragment_mz > 67.048 & fragment_mz < 67.06]
  cat(sprintf("p0=%2d: %d group(s) near 67.054 Da (samples: %s)\n",
              p0, nrow(near67), paste(near67$n_samples, collapse = ", ")))
}
#> p0= 4: 2 group(s) near 67.054 Da (samples: 5, 6)
#> p0= 8: 1 group(s) near 67.054 Da (samples: 10)
#> p0=12: 1 group(s) near 67.054 Da (samples: 10)
```

At `p0 = 4`, the fragment near 67.054 Da is reported as two separate
groups ~10 ppm apart, each seen in only part of the samples:
between-sample calibration drift, wider than a `p0` sized for within-run
precision alone. From `p0 = 8` they merge into one group seen in all 10
samples. This is exactly the kind of thing to check when tuning a
grouper for a new instrument or acquisition mode: pick a handful of
real, multi-sample features and look for this pattern before trusting a
`p0` value across a whole cohort.

## Writing your own grouper

A grouper is a factory function returning a closure with a fixed
signature — see `?grouper-contract`.
\[[`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md)\]
calls it once per feature with `mz`, `intensity` and `target_mz`.

``` r

# A trivial custom grouper: fixed ppm tolerance, but reject a merge if the
# combined evidence is too thin (n < 3) - contrived, just to show the shape.
sparse_aware_grouper <- function(ppm = 20, min_n = 3L) {
  force(ppm); force(min_n)
  function(mz, intensity = NULL, target_mz = NULL, scan = NULL, ms1_int = NULL) {
    if (length(mz) < min_n) return(seq_along(mz))  # too little evidence: don't merge anything
    fast_1D_upgma_grouping_cpp(mz, maxDiff = 1e-6 * target_mz * ppm)
  }
}

validate_grouper(sparse_aware_grouper())

stats_custom <- process_pseudo_spec(raw, mz_target = mz_target,
                                    grouper = sparse_aware_grouper(ppm = 20))
nrow(stats_custom)
#> [1] 633
```

\[[`validate_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/validate_grouper.md)\]
sanity-checks the contract (right output length, no more groups than
distinct input values) against small synthetic input before you run it
against real data — do this once when wiring up something new.

## Writing your own filter

Statistics (`group_statistics_fn`) and filtering (`filter_fn`) are
separate closures for the same reason grouping is swappable — see
`?stats-contract` and `?filter-contract`. The default pair is
\[[`correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_statistics.md)\] +
\[[`correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_filtering.md)\];
\[[`apex_correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_statistics.md)\] +
\[[`apex_correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_filtering.md)\]
add intra-sample MS1↔︎MS2 apex-agreement columns and thresholds, for when
correlation alone lets through fragments from a coeluting compound.

``` r

# A stricter filter: require both a higher correlation AND more supporting
# (sample, ms1_scan) cells than the package default.
strict_filter <- function(stats) {
  stats[correlation > 0.8 & count >= 30]
}

clean_pseudo_spectra(
  targets        = "FT00964",
  file_directory = fixture_dir,
  fd             = fd,
  grouper        = iaac_grouper(k = 2, p0 = 12, Isat = 1e4),
  filter_fn      = strict_filter,
  report         = FALSE,
  mgf_dir        = NULL,
  progress       = FALSE
)[, .(feature_id, mz_g2, count, r, fdr)]
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpehZATr/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#>    feature_id     mz_g2 count         r          fdr
#>        <char>     <num> <int>     <num>        <num>
#> 1:    FT00964  67.05413    96 0.9168420 1.463360e-36
#> 2:    FT00964 138.05473    96 0.8754481 4.833028e-29
```

## Where this is checked automatically

`tests/testthat/test-process_pseudo_spec.R` pins the numbers shown above
(the `p0 = 12` top-group correlation, the `p0 = 4` vs `p0 = 12`
split/merge behavior) as regression tests against this same fixture. If
you change the grouping design and those tests start failing, that’s
either a real regression or a deliberate change worth calling out
explicitly — not something to silently update without checking why the
numbers moved.
