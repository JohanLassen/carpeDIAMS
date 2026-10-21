# Build a pseudo-spectra MGF (dispatcher)

Picks
[`process_xcms()`](https://johanlassen.github.io/carpeDIAMS/reference/process_xcms.md),
[`process_targetlist()`](https://johanlassen.github.io/carpeDIAMS/reference/process_targetlist.md)
or
[`process_mzml()`](https://johanlassen.github.io/carpeDIAMS/reference/process_mzml.md)
from `input` and `targetlist`, and passes everything else on.

## Usage

``` r
carpediams(
  input,
  output,
  targetlist = NULL,
  tmp_dir = NULL,
  grouper = NULL,
  group_statistics_fn = NULL,
  filter_fn = NULL,
  ...
)
```

## Arguments

- input:

  Object or path; see Routing.

- output:

  MGF path. Not needed for the `"prepare"`, `"sample"` and `"filter"`
  stages.

- targetlist:

  `NULL`, or a CSV path / data.frame with `mz` and `rt` columns.

- tmp_dir:

  Working directory. `NULL` uses a
  [`tempfile()`](https://rdrr.io/r/base/tempfile.html) that is **deleted
  when the call returns**, so staged runs must pass a path.

- grouper:

  A grouper closure; `NULL` uses the route default
  (`ppm_upgma(ppm = 20)` for xcms/target-list via their `ppm` argument,
  `ppm_upgma(ppm = 30)` for mzML).

- group_statistics_fn:

  A statistics closure; `NULL` uses
  `correlation_statistics(min_count = 10)`.

- filter_fn:

  A filter closure; `NULL` uses
  `correlation_filtering(cor_thresh = 0.1, fdr_thresh = 0.1)`.

- ...:

  Passed to the routed function (`stage`, `ppm`, `top_n`, …).
  [`process_mzml()`](https://johanlassen.github.io/carpeDIAMS/reference/process_mzml.md)
  takes no extra arguments.

## Value

What the routed function returns: `tmp_dir` (invisibly) for the xcms and
target-list routes, the pseudo-spectra table for mzML.

## Routing

|  |  |  |
|----|----|----|
| `input` | `targetlist` | Route |
| `XcmsExperiment`/`XCMSnExp`, or `.rds`/`.qs`/`.rda` path | ignored (warning) | [`process_xcms()`](https://johanlassen.github.io/carpeDIAMS/reference/process_xcms.md) |
| mzML file, mzML directory, or vector of paths | given | [`process_targetlist()`](https://johanlassen.github.io/carpeDIAMS/reference/process_targetlist.md) |
| one mzML file | `NULL` | [`process_mzml()`](https://johanlassen.github.io/carpeDIAMS/reference/process_mzml.md) |
| mzML directory or vector of paths | `NULL` | error |

## Pluggable steps

- `grouper`
  ([grouper-contract](https://johanlassen.github.io/carpeDIAMS/reference/grouper-contract.md)):
  [`ppm_upgma()`](https://johanlassen.github.io/carpeDIAMS/reference/ppm_upgma.md),
  [`abs_upgma()`](https://johanlassen.github.io/carpeDIAMS/reference/abs_upgma.md),
  [`iaac_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/iaac_grouper.md),
  [`pelt_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/pelt_grouper.md),
  [`bin_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/bin_grouper.md).

- `group_statistics_fn`
  ([stats-contract](https://johanlassen.github.io/carpeDIAMS/reference/stats-contract.md)):
  [`correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_statistics.md),
  [`apex_correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_statistics.md).

- `filter_fn`
  ([filter-contract](https://johanlassen.github.io/carpeDIAMS/reference/filter-contract.md)):
  [`correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_filtering.md),
  [`apex_correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_filtering.md)
  (needs
  [`apex_correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_statistics.md)).

- `consensus_mz_fn`
  ([consensus-mz-contract](https://johanlassen.github.io/carpeDIAMS/reference/consensus-mz-contract.md)):
  [`consensus_mz_weighted()`](https://johanlassen.github.io/carpeDIAMS/reference/consensus_mz_weighted.md).
  Pass it through `...`; xcms and target-list routes only.

## See also

Other api:
[`process_mzml()`](https://johanlassen.github.io/carpeDIAMS/reference/process_mzml.md),
[`process_targetlist()`](https://johanlassen.github.io/carpeDIAMS/reference/process_targetlist.md),
[`process_xcms()`](https://johanlassen.github.io/carpeDIAMS/reference/process_xcms.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# default — correlation-based filter
carpediams(input = "grouped.rds", output = "out.mgf")

# swap the grouper
carpediams(input = "grouped.rds", output = "out.mgf",
           grouper = iaac_grouper(k = 1, p0 = 8, Isat = 1e5))

# apex-agreement filter (needs the matching statistics)
carpediams(input = "grouped.rds", output = "out.mgf",
           group_statistics_fn = apex_correlation_statistics(),
           filter_fn           = apex_correlation_filtering(
                                   max_apex_shift_med = 0.5,
                                   max_apex_shift_iqr = 0.75))
} # }
```
