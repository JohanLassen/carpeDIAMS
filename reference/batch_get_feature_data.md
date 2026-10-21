# Fetch raw MS2 + MS1 rows for a vector of features from DuckDB

Pairs each MS1 pointer row with the MS2 centroids of the **next scan**
(`pool.scan = pointers.scan + 1`, same sample), keeping fragments with
m/z at most precursor + 2. MS1 rows with no matching MS2 centroid are
kept with `NA` fragment columns.

## Usage

``` r
batch_get_feature_data(feature_ids, con, chunk_size = 500)
```

## Arguments

- feature_ids:

  Character vector of feature ids.

- con:

  A connection from
  [`setup_connection()`](https://johanlassen.github.io/carpeDIAMS/reference/setup_connection.md).

- chunk_size:

  Features per query.

## Value

`data.table`, one row per (MS1 scan, MS2 centroid):
`sample, feature, mz` (precursor),
`fragment_mz, int, ms1_scan, ms2_scan, ms1_int, apex_scan, ms1_rt, ms2_rt`.
`apex_scan` is the apex found upstream (xcms peak or nearest scan to the
target RT); the RT columns are `NULL`/`NA` for parquets written without
`rt`.
