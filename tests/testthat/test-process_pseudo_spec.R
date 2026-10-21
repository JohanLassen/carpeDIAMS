# Regression tests for process_pseudo_spec()'s grouping, built on a real,
# previously-problematic feature rather than synthetic data.
#
# Fixture: inst/testdata/FT00964/ - a faithful (not synthetic) slice of a
# real carpeDIAMS production run's pool/pointers parquets, trimmed to
# exactly the rows the real SQL join in batch_get_feature_data() ever
# touches for feature "FT00964" (10 samples).
#
# NEEDS RE-PIN: every count/correlation/fragment_mz value below was measured
# under the old two-pass (within-sample cluster -> broadcast consensus ->
# cross-sample cluster) design. process_pseudo_spec() now runs a single
# pooled-clustering pass with no broadcast/denoise step - a different enough
# algorithm that these numbers are not expected to still hold. Re-run on the
# cluster (this sandbox has no duckdb/arrow) and re-pin.

fixture_dir <- system.file("testdata", "FT00964", package = "carpeDIAMS")
fd_fixture  <- readRDS(file.path(fixture_dir, "fd.rds"))

test_that("IAAC grouping (production p0=12) reproduces the known-good real result", {
  testthat::skip_if_not_installed("duckdb")
  testthat::skip_if_not_installed("arrow")

  con <- setup_connection(fixture_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  dt <- batch_get_feature_data("FT00964", con)
  expect_gt(nrow(dt), 0)
  expect_setequal(unique(dt$sample), c(6, 7, 8, 11, 14, 18, 19, 24, 25, 27))

  mz_target <- fd_fixture["FT00964", "mzmed"]

  # This is the actual production grouper config for the xcms route
  # (08_carpeDIAMS/drugs_workflow/scripts/01/02/03_*_xcms.R): the route
  # this exact p0 was already validated against at plan time.
  grouper <- iaac_grouper(k = 2, p0 = 12, Isat = 1e4)
  stats <- process_pseudo_spec(dt, mz_target = mz_target, grouper = grouper)

  expect_true(nrow(stats) > 0)
  top <- stats[order(-abs(correlation))][1]

  # Real, previously-verified numbers (see the fixture header comment).
  expect_equal(top$count, 96L)
  expect_equal(top$correlation, 0.9168420, tolerance = 1e-4)
  expect_equal(top$fragment_mz, 67.05413, tolerance = 1e-3)
})

test_that("IAAC grouping does not under-merge at a p0 wide enough for this feature's between-sample drift", {
  testthat::skip_if_not_installed("duckdb")
  testthat::skip_if_not_installed("arrow")

  con <- setup_connection(fixture_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  dt <- batch_get_feature_data("FT00964", con)
  mz_target <- fd_fixture["FT00964", "mzmed"]

  # At p0=4 (the target-list route's old production value, since retuned -
  # see 00_process_targetlist.R), a real fragment near 67.054 Da splits
  # into two sub-clusters ~10 ppm apart instead of merging - the concrete
  # example that motivated retuning that route's p0. This test pins the
  # boundary: it should NOT be merged at p0=4, and SHOULD be merged (one
  # group, near-full sample support) at p0=12 (the now-shared value).
  grouper_tight <- iaac_grouper(k = 1, p0 = 4, Isat = 1e5)
  stats_tight <- process_pseudo_spec(dt, mz_target = mz_target, grouper = grouper_tight)
  near67_tight <- stats_tight[fragment_mz > 67.04 & fragment_mz < 67.06]
  expect_gt(nrow(near67_tight), 1)  # still split at p0=4

  grouper_wide <- iaac_grouper(k = 1, p0 = 12, Isat = 1e5)
  stats_wide <- process_pseudo_spec(dt, mz_target = mz_target, grouper = grouper_wide)
  near67_wide <- stats_wide[fragment_mz > 67.048 & fragment_mz < 67.06]
  # merged into (at most) one well-supported group near 67.054
  expect_lte(nrow(near67_wide[count > 50]), 1)
  if (nrow(near67_wide[count > 50]) == 1) {
    expect_gte(near67_wide[count > 50]$count, 80L)
  }
})

test_that("process_pseudo_spec() returns a well-formed per-group table with no NA group labels", {
  testthat::skip_if_not_installed("duckdb")
  testthat::skip_if_not_installed("arrow")

  con <- setup_connection(fixture_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  dt <- batch_get_feature_data("FT00964", con)
  mz_target <- fd_fixture["FT00964", "mzmed"]

  stats <- process_pseudo_spec(dt, mz_target = mz_target,
                               grouper = ppm_upgma(ppm = 30))
  expect_true(nrow(stats) > 0)
  expect_false(anyNA(stats$group))
  expect_true(all(c("fragment_mz", "correlation", "count", "n_samples") %in% colnames(stats)))
})

test_that("clean_pseudo_spectra() end-to-end smoke test on the real fixture", {
  testthat::skip_if_not_installed("duckdb")
  testthat::skip_if_not_installed("arrow")

  cleaned <- clean_pseudo_spectra(
    targets        = "FT00964",
    file_directory = fixture_dir,
    fd             = fd_fixture,
    grouper        = iaac_grouper(k = 2, p0 = 12, Isat = 1e4),
    report         = FALSE,
    mgf_dir        = NULL,
    progress       = FALSE
  )

  expect_true(nrow(cleaned) > 0)
  expect_true(all(cleaned$feature_id == "FT00964"))
  # deconv_error is NA on every row - this feature does produce real spectra
  expect_true(all(is.na(cleaned$deconv_error)))

  top <- cleaned[order(-abs(r))][1]
  expect_equal(top$count, 96L)
  expect_equal(top$r, 0.9168420, tolerance = 1e-3)
})

test_that("apex_source/window_mode defaults are bit-identical to the pre-existing behavior", {
  testthat::skip_if_not_installed("duckdb")
  testthat::skip_if_not_installed("arrow")

  con <- setup_connection(fixture_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  dt <- batch_get_feature_data("FT00964", con)
  mz_target <- fd_fixture["FT00964", "mzmed"]
  grouper   <- iaac_grouper(k = 2, p0 = 12, Isat = 1e4)

  stats_implicit <- process_pseudo_spec(dt, mz_target = mz_target, grouper = grouper)
  stats_explicit <- process_pseudo_spec(dt, mz_target = mz_target, grouper = grouper,
                                        apex_source = "recompute", window_mode = "buffer")
  expect_equal(stats_implicit, stats_explicit)

  top <- stats_implicit[order(-abs(correlation))][1]
  expect_equal(top$count, 96L)
  expect_equal(top$correlation, 0.9168420, tolerance = 1e-4)
  expect_equal(top$fragment_mz, 67.05413, tolerance = 1e-3)
})

test_that("batch_get_feature_data() carries apex_scan, constant per (sample, feature)", {
  testthat::skip_if_not_installed("duckdb")
  testthat::skip_if_not_installed("arrow")

  con <- setup_connection(fixture_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  dt <- batch_get_feature_data("FT00964", con)
  expect_true("apex_scan" %in% colnames(dt))
  expect_true(dt[, all(data.table::uniqueN(apex_scan) == 1), by = sample]$V1 |> all())

  # Real, previously-verified divergence on this fixture (see the 2026-09
  # apex-source investigation): several samples' upstream apex_scan differs
  # from ms1_scan[which.max(ms1_int)] - confirms the two apex_source modes
  # aren't a no-op on real data.
  local_apex    <- dt[!is.na(int) & !is.na(fragment_mz), .(local = ms1_scan[which.max(ms1_int)]), by = sample]
  upstream_apex <- dt[, .(upstream = apex_scan[1]), by = sample]
  cmp <- merge(local_apex, upstream_apex, by = "sample")
  expect_gt(sum(cmp$local != cmp$upstream), 0)
})

test_that("apex_source changes the result under window_mode = 'fwhm' (where apex precision matters)", {
  testthat::skip_if_not_installed("duckdb")
  testthat::skip_if_not_installed("arrow")

  con <- setup_connection(fixture_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  dt <- batch_get_feature_data("FT00964", con)
  mz_target <- fd_fixture["FT00964", "mzmed"]
  grouper   <- iaac_grouper(k = 2, p0 = 12, Isat = 1e4)

  # Under the default wide buffer (50 scans), a few scans' worth of apex
  # disagreement doesn't change which rows fall in the window - apex_source
  # is a genuine no-op here on this fixture. FWHM's much tighter window is
  # where the choice of apex actually matters; this pins that real
  # difference rather than a no-op.
  stats_recompute <- process_pseudo_spec(dt, mz_target = mz_target, grouper = grouper,
                                         window_mode = "fwhm", apex_source = "recompute")
  stats_upstream  <- process_pseudo_spec(dt, mz_target = mz_target, grouper = grouper,
                                         window_mode = "fwhm", apex_source = "upstream")
  expect_false(isTRUE(all.equal(stats_recompute[order(fragment_mz)], stats_upstream[order(fragment_mz)])))

  top_recompute <- stats_recompute[fdr < 0.05][order(-abs(correlation))][1]
  expect_equal(top_recompute$count, 30L)
  expect_equal(top_recompute$correlation, 0.6723145, tolerance = 1e-4)
  expect_equal(top_recompute$fragment_mz, 67.05409, tolerance = 1e-3)

  top_upstream <- stats_upstream[order(-abs(correlation))][1]
  expect_equal(top_upstream$count, 11L)
  expect_equal(top_upstream$correlation, 0.765987, tolerance = 1e-4)
  expect_equal(top_upstream$fragment_mz, 154.9816, tolerance = 1e-3)
})

test_that(".fwhm_bounds() is anchored at the apex's own intensity, not the trace's max", {
  # Clean symmetric peak: apex at scan 6 (int 120), half-max = 60.
  scan <- 1:11
  int  <- c(10, 20, 40, 70, 100, 120, 100, 70, 40, 20, 10)
  b <- carpeDIAMS:::.fwhm_bounds(scan, int, apex_scan = 6)
  expect_equal(b$lo, 4); expect_equal(b$hi, 8)

  # A co-eluting neighbor far away (scan 15) has intensity 200 - HIGHER than
  # the true apex (scan 6, int 50). The threshold must stay anchored at the
  # apex's own intensity (half = 25), not the trace's max (half = 100) - the
  # whole point of the mechanism. The neighbor must not appear in the window.
  scan2 <- c(1:11, 15)
  int2  <- c(5, 10, 20, 35, 50, 50, 35, 20, 10, 5, 2, 200)
  b2 <- carpeDIAMS:::.fwhm_bounds(scan2, int2, apex_scan = 6)
  expect_equal(b2$lo, 4); expect_equal(b2$hi, 7)
  expect_true(b2$hi < 15)

  # Monotonically increasing trace that never dips below half-max of its
  # last point within the available data - clamp to the data's edge, don't
  # error, extrapolate, or return an unbounded window.
  scan3 <- 1:5
  int3  <- c(100, 105, 110, 115, 120)
  b3 <- carpeDIAMS:::.fwhm_bounds(scan3, int3, apex_scan = 5)
  expect_equal(b3$lo, 1); expect_equal(b3$hi, 5)
})

test_that("clean_pseudo_spectra() runs end-to-end with both new params non-default", {
  testthat::skip_if_not_installed("duckdb")
  testthat::skip_if_not_installed("arrow")

  cleaned <- clean_pseudo_spectra(
    targets        = "FT00964",
    file_directory = fixture_dir,
    fd             = fd_fixture,
    grouper        = iaac_grouper(k = 2, p0 = 12, Isat = 1e4),
    apex_source    = "upstream",
    window_mode    = "fwhm",
    report         = FALSE,
    mgf_dir        = NULL,
    progress       = FALSE
  )

  expect_true(nrow(cleaned) > 0)
  expect_true(all(cleaned$feature_id == "FT00964"))
})

test_that("consensus_mz_fn = NULL is bit-identical to the pre-existing hardcoded default", {
  testthat::skip_if_not_installed("duckdb")
  testthat::skip_if_not_installed("arrow")

  con <- setup_connection(fixture_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  dt <- batch_get_feature_data("FT00964", con)
  mz_target <- fd_fixture["FT00964", "mzmed"]
  grouper   <- iaac_grouper(k = 2, p0 = 12, Isat = 1e4)

  # consensus_mz_fn = NULL must resolve to exactly today's hardcoded
  # frac=0.5, Isat=1e4 behavior - same pinned numbers as the very first
  # test in this file, now also reachable via an explicit closure.
  stats_default  <- process_pseudo_spec(dt, mz_target = mz_target, grouper = grouper)
  stats_explicit <- process_pseudo_spec(dt, mz_target = mz_target, grouper = grouper,
                                        consensus_mz_fn = consensus_mz_weighted(frac = 0.5, Isat = 1e4))

  expect_equal(stats_default, stats_explicit)

  top <- stats_default[order(-abs(correlation))][1]
  expect_equal(top$count, 96L)
  expect_equal(top$correlation, 0.9168420, tolerance = 1e-4)
  expect_equal(top$fragment_mz, 67.05413, tolerance = 1e-3)
})
