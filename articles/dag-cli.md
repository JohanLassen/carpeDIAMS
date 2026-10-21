# Running stages as cluster jobs

The xcms and target-list routes can be split into four stages that only
communicate through files in a shared `tmp_dir`. Each stage is one
[`carpediams()`](https://johanlassen.github.io/carpeDIAMS/reference/carpediams.md)
call, so a workflow manager (gwf, Snakemake, Nextflow, SLURM arrays) can
run the per-sample and per-batch stages in parallel.

    prepare ──> sample 1..N ──> filter 1..B ──> finalize ──> MGF

For a step-by-step version in one R session, see
[`vignette("xcms-interactive")`](https://johanlassen.github.io/carpeDIAMS/articles/xcms-interactive.md).

## Stages

| Stage | Runs | Reads | Writes |
|:---|:---|:---|:---|
| `prepare` | once | `input` (and `targetlist`) | `manifest.csv`, `fd.rds`; xcms route also `mssubset.rds`, `sample_topfeatures.rds` |
| `sample` | per file | files from `prepare`, the mzML | `sample_<i>_pool.parquet`, `sample_<i>_pointers.parquet` |
| `filter` | per batch | `fd.rds`, all sample parquets | `features_batch_<b>.parquet` |
| `finalize` | once | `fd.rds`, all `features_batch_*.parquet` | the MGF |

`<i>` is `sample_idx` in `manifest.csv`. `filter` splits the features
into `n_batches` contiguous slices and processes slice `batch_idx`.
`finalize` binds every `features_batch_*.parquet` it finds. If there are
none, it filters all features itself, which makes the `filter` stage
optional for small runs.

## Arguments every call needs

- `input`:
  [`carpediams()`](https://johanlassen.github.io/carpeDIAMS/reference/carpediams.md)
  uses it to pick the route. Only `prepare` reads it, so pass the path,
  not a loaded object.
- `targetlist` (target-list route): needed at **every** stage. Without
  it, a single mzML input is routed to
  [`process_mzml()`](https://johanlassen.github.io/carpeDIAMS/reference/process_mzml.md).
- `tmp_dir`: the same directory for every stage. If it is left as
  `NULL`, a temporary directory is created and deleted when the call
  returns.
- `output`: only needed by `finalize`.

Keep stage-specific parameters identical across jobs: `ppm` and
`rt_buffer` (xcms) or `ppm` and `rt_width` (target list) for every
`sample` job, and the same `grouper`, `group_statistics_fn` and
`filter_fn` for every `filter` job and for `finalize`.

## The calls (xcms route)

``` r

library(carpeDIAMS)
input   <- "grouped_xcms.rds"
tmp_dir <- "work/carpediams"

# prepare
carpediams(input = input, tmp_dir = tmp_dir, stage = "prepare",
           subsample = TRUE, subset_size = 100, top_n = 10)

# sample, one job per row of manifest.csv
carpediams(input = input, tmp_dir = tmp_dir, stage = "sample",
           sample_idx = 1, ppm = 20, rt_buffer = 5)

# filter, one job per batch
carpediams(input = input, tmp_dir = tmp_dir, stage = "filter",
           batch_idx = 1, n_batches = 200,
           grouper   = pelt_grouper(p0 = 15, Isat = 1e5),
           filter_fn = correlation_filtering(cor_thresh = 0.5, fdr_thresh = 0.01))

# finalize
carpediams(input = input, tmp_dir = tmp_dir, stage = "finalize",
           output = "results/all_features.mgf",
           grouper   = pelt_grouper(p0 = 15, Isat = 1e5),
           filter_fn = correlation_filtering(cor_thresh = 0.5, fdr_thresh = 0.01))
```

`sample_file = "<basename>.mzML"` can be used in place of `sample_idx`.

The target-list route is the same, with `targetlist` added to every call
and `rt_width` in place of `rt_buffer`. Its `prepare` has no `subsample`
or `top_n`, because every target is extracted from every file.

## From the shell

Each call fits in an `Rscript -e` line, which is all a job template
needs:

``` sh
Rscript -e 'carpeDIAMS::carpediams(input = "grouped_xcms.rds",
  tmp_dir = "work/carpediams", stage = "sample", sample_idx = 17,
  ppm = 20, rt_buffer = 5)'
```

The number of `sample` jobs is the number of rows in `manifest.csv`,
which only exists after `prepare`. Either read it in a second planning
phase (Snakemake checkpoint), or declare the jobs up front from a known
sample count. `n_batches` can be any number, because batches past the
last feature write an empty parquet.

## Memory

- `sample`: holds one file’s spectra.
- `filter`: DuckDB reads the parquets of all samples for the batch’s
  features. It is capped at `memory_limit = '128GB'` and 2 threads, and
  spills to `<tmp_dir>/.duckdb_tmp`. Raise `n_batches` to lower memory
  per job.
- `finalize` with batches: reads the batch parquets, which are small.

## Failures and reruns

- A `sample` or `filter` job failed: rerun that job; it overwrites its
  own files.
- Changed filter settings: delete `features_batch_*.parquet`, rerun
  `filter` and `finalize`. The sample parquets are reused.
- Changed `n_batches`: delete the old batch files first, because
  `finalize` binds every batch file it finds.
- Changed `prepare` or `sample` settings: start from an empty `tmp_dir`.

Checking that all sample files exist before filtering:

``` r

manifest <- read.csv(file.path(tmp_dir, "manifest.csv"))
expected <- c(sprintf("sample_%d_pool.parquet",     manifest$sample_idx),
              sprintf("sample_%d_pointers.parquet", manifest$sample_idx))
setdiff(expected, list.files(tmp_dir))
```

On the xcms route, a sample with no top-`top_n` feature writes no files,
so a missing pair is not always an error. Look for the sample’s message
in the job log.
