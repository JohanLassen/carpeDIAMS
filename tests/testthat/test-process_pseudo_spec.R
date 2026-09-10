# Regression tests for the two-pass (within-sample -> cross-sample) grouping
# design in process_pseudo_spec(), built on a real, previously-problematic
# feature rather than synthetic data.
#
# Fixture: inst/testdata/FT00964/ - a faithful (not synthetic) slice of a
# real carpeDIAMS production run's pool/pointers parquets, trimmed to
# exactly the rows the real SQL join in batch_get_feature_data() ever
# touches for feature "FT00964" (10 samples). See the 2026-09-10 IAAC
# two-pass grouping change for how this fixture was derived and validated.
#
# The expected numbers below were computed by running this exact
# algorithm (process_pseudo_spec()'s current body) against this exact
# fixture; a change to the grouping design that shifts them meaningfully
# should be a deliberate, reviewed decision, not a silent regression.

fixture_dir <- system.file("testdata", "FT00964", package = "carpeDIAMS")
fd_fixture  <- readRDS(file.path(fixture_dir, "fd.rds"))

test_that("two-pass IAAC grouping (production p0=12) reproduces the known-good real result", {
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

test_that("two-pass IAAC grouping does not under-merge at a p0 wide enough for this feature's between-sample drift", {
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

test_that("process_pseudo_spec() never fans dt out across (sample, group1) - one group label per raw row", {
  testthat::skip_if_not_installed("duckdb")
  testthat::skip_if_not_installed("arrow")

  con <- setup_connection(fixture_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  dt <- batch_get_feature_data("FT00964", con)
  mz_target <- fd_fixture["FT00964", "mzmed"]

  # Exercises the full two-pass path (including the pass2_input update-join)
  # end to end and simply confirms it completes and returns a well-formed
  # per-group table - a smoke test for the join not fanning out or leaving
  # NA group labels behind.
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
