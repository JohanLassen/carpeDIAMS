# Package index

## Build an MGF

Entry points.
[`carpediams()`](https://johanlassen.github.io/carpeDIAMS/reference/carpediams.md)
picks one of the three routes.

- [`carpediams()`](https://johanlassen.github.io/carpeDIAMS/reference/carpediams.md)
  : Build a pseudo-spectra MGF (dispatcher)
- [`process_xcms()`](https://johanlassen.github.io/carpeDIAMS/reference/process_xcms.md)
  : Build a pseudo-spectra MGF from a grouped xcms object
- [`process_targetlist()`](https://johanlassen.github.io/carpeDIAMS/reference/process_targetlist.md)
  : Build a pseudo-spectra MGF from mzML file(s) and a target list
- [`process_mzml()`](https://johanlassen.github.io/carpeDIAMS/reference/process_mzml.md)
  : Build a pseudo-spectra MGF from one mzML file

## Fragment grouping

Grouper factories. Pass one as `grouper`.

- [`grouper-contract`](https://johanlassen.github.io/carpeDIAMS/reference/grouper-contract.md)
  : Grouper contract
- [`ppm_upgma()`](https://johanlassen.github.io/carpeDIAMS/reference/ppm_upgma.md)
  : ppm-based UPGMA grouper (the package default)
- [`abs_upgma()`](https://johanlassen.github.io/carpeDIAMS/reference/abs_upgma.md)
  : Absolute-tolerance UPGMA grouper
- [`iaac_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/iaac_grouper.md)
  : Intensity-adaptive agglomerative clustering grouper
- [`pelt_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/pelt_grouper.md)
  : PELT (penalized optimal partitioning) grouper
- [`bin_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/bin_grouper.md)
  : Fixed-width binning grouper
- [`validate_grouper()`](https://johanlassen.github.io/carpeDIAMS/reference/validate_grouper.md)
  : Sanity-check a custom grouper against the contract

## Statistics and filters

Score fragment groups (`group_statistics_fn`) and keep the good ones
(`filter_fn`).

- [`stats-contract`](https://johanlassen.github.io/carpeDIAMS/reference/stats-contract.md)
  : Statistics contract
- [`correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_statistics.md)
  : Per-group correlation statistics (the package default)
- [`apex_correlation_statistics()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_statistics.md)
  : Per-group correlation + apex-agreement statistics
- [`filter-contract`](https://johanlassen.github.io/carpeDIAMS/reference/filter-contract.md)
  : Filter contract
- [`correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/correlation_filtering.md)
  : Correlation + BH FDR filter (the package default)
- [`apex_correlation_filtering()`](https://johanlassen.github.io/carpeDIAMS/reference/apex_correlation_filtering.md)
  : Correlation + FDR + apex-agreement filter

## Consensus m/z

- [`consensus-mz-contract`](https://johanlassen.github.io/carpeDIAMS/reference/consensus-mz-contract.md)
  : Consensus m/z contract
- [`consensus_mz_weighted()`](https://johanlassen.github.io/carpeDIAMS/reference/consensus_mz_weighted.md)
  : Intensity-weighted FWHM-core consensus m/z (the package default)

## Spectrum building

The per-feature machinery behind the entry points.

- [`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md)
  : Clean pseudo-spectra for many features (main batch entry)
- [`process_pseudo_spec()`](https://johanlassen.github.io/carpeDIAMS/reference/process_pseudo_spec.md)
  : Build a feature's per-group statistics table from raw MS1+MS2
  data.table
- [`export_mgf()`](https://johanlassen.github.io/carpeDIAMS/reference/export_mgf.md)
  : Assemble a combined MGF from a flat pseudo-spectra results table
- [`export_pseudospectra()`](https://johanlassen.github.io/carpeDIAMS/reference/export_pseudospectra.md)
  : Write a list of spectra as one MGF and optional SIRIUS .ms files

## Low-level

- [`fast_1D_upgma_grouping_cpp()`](https://johanlassen.github.io/carpeDIAMS/reference/fast_1D_upgma_grouping_cpp.md)
  : Fast 1D UPGMA grouping
- [`fast_1D_iaac_grouping_cpp()`](https://johanlassen.github.io/carpeDIAMS/reference/fast_1D_iaac_grouping_cpp.md)
  : Fast 1D intensity-adaptive grouping (IAAC), O(n log n)
- [`fast_1D_pelt_grouping_cpp()`](https://johanlassen.github.io/carpeDIAMS/reference/fast_1D_pelt_grouping_cpp.md)
  : Exact 1D grouping via PELT (penalized optimal partitioning)
