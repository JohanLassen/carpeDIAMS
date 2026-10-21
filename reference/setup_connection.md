# Open a DuckDB connection with views over pool/pointer parquets

Opens an in-memory DuckDB with views `pointers` and `pool` over all
`*_pointers.parquet` / `*_pool.parquet` files in `file_directory`.
DuckDB spills to `<file_directory>/.duckdb_tmp`, with a fixed
`memory_limit = '128GB'` and `threads = 2`. Errors if there are no
pointer files. Close it with `DBI::dbDisconnect(con, shutdown = TRUE)`;
[`clean_pseudo_spectra()`](https://johanlassen.github.io/carpeDIAMS/reference/clean_pseudo_spectra.md)
does not close a connection you pass in.

## Usage

``` r
setup_connection(file_directory)
```

## Arguments

- file_directory:

  Directory with the per-sample parquets.

## Value

A DuckDB connection.
