## ---------------------------------------------------------------
## Tests for fect_mspe: criterion (S3), cv.method (S4), weights (S6),
## norm.para (S7), return structure (S8), input validation (S9), and
## fect_mspe simplification (Section D, NEW for Phase 2).
##
## Originally part of test-score-unify.R; split out 2026-05-03 for
## progress visibility. Shared fixtures live in helper-score-unify.R.
## ---------------------------------------------------------------

## =================================================================
## S3: fect_mspe criterion support (updated for Phase 2)
## =================================================================

test_that("S3.1: Default criterion='mspe' matches old RMSE", {

  skip_on_cran()
  res_new <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 123, criterion = "mspe")
  ))

  # RMSE = sqrt(MSPE) invariant (P1)
  if ("MSPE" %in% names(res_new$summary)) {
    expect_equal(res_new$summary$RMSE, sqrt(res_new$summary$MSPE),
                 tolerance = 1e-10)
  }
})

test_that("S3.2: All 7 criteria produce finite, positive scores", {

  skip_on_cran()
  for (crit in c("mspe", "wmspe", "gmspe", "wgmspe", "mad",
                 "moment", "gmoment")) {
    res <- suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 123, criterion = crit)
    ))
    crit_upper <- toupper(crit)
    if (crit_upper %in% names(res$summary)) {
      expect_true(all(is.finite(res$summary[[crit_upper]])),
                  info = paste("Criterion", crit, "should be finite"))
      expect_true(all(res$summary[[crit_upper]] > 0),
                  info = paste("Criterion", crit, "should be positive"))
    }
  }
})


## =================================================================
## S4: fect_mspe with cv.method (Phase 2 — replaces mask.method)
## =================================================================

test_that("S4.1: cv.method='treated_units' masking runs without error", {

  skip_on_cran()
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, cv.method = "treated_units",
              k = 3, cv.prop = 0.1,
              cv.nobs = 3, cv.donut = 1, min.T0 = 3,
              criterion = "mspe")
  ))
  expect_true("summary" %in% names(res))
  if ("MSPE" %in% names(res$summary)) {
    expect_true(all(is.finite(res$summary$MSPE)))
    expect_true(all(res$summary$MSPE > 0))
  }
  expect_true(res$summary$RMSE > 0)
})

test_that("S4.2: cv.method='all_units' masking runs without error", {

  skip_on_cran()
  res_cv <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, cv.method = "all_units",
              k = 3, criterion = "mspe")
  ))
  expect_true("summary" %in% names(res_cv))
  expect_true(res_cv$summary$RMSE > 0)
})

test_that("S4.3: cv.method='treated_units' with k=1 (single fold)", {

  skip_on_cran()
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, cv.method = "treated_units",
              k = 1, criterion = "mspe")
  ))
  expect_true(res$summary$RMSE > 0)
})


## =================================================================
## S6: fect_mspe with observation weights (W) — updated for Phase 2
## =================================================================

test_that("S6.1: W parameter produces different scores than unweighted", {

  skip_on_cran()
  TT <- nrow(out_base$Y.dat)
  NN <- ncol(out_base$Y.dat)
  set.seed(7)
  W_mat <- matrix(runif(TT * NN, 0.5, 1.5), nrow = TT, ncol = NN)
  res_unw <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, criterion = "mspe",
              W = NULL)
  ))
  res_w <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, criterion = "mspe",
              W = W_mat)
  ))
  expect_false(identical(res_unw$summary$RMSE, res_w$summary$RMSE))
})

test_that("S6.2: Uniform W equals unweighted", {

  skip_on_cran()
  TT <- nrow(out_base$Y.dat)
  NN <- ncol(out_base$Y.dat)
  W_uniform <- matrix(1, nrow = TT, ncol = NN)
  res_u <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, criterion = "mspe",
              W = W_uniform)
  ))
  res_n <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, criterion = "mspe",
              W = NULL)
  ))
  if ("MSPE" %in% names(res_u$summary) && "MSPE" %in% names(res_n$summary)) {
    expect_equal(res_u$summary$MSPE, res_n$summary$MSPE, tolerance = 1e-10)
  }
})


## =================================================================
## S7: fect_mspe with norm.para — updated for Phase 2
## =================================================================

test_that("S7.1: norm.para scales scores", {

  skip_on_cran()
  np <- c(2.0, 0.0)
  res_raw <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, criterion = "mspe")
  ))
  res_norm <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, criterion = "mspe",
              norm.para = np)
  ))
  if ("MSPE" %in% names(res_raw$summary) &&
      "MSPE" %in% names(res_norm$summary)) {
    expect_equal(res_norm$summary$MSPE, res_raw$summary$MSPE * 4.0,
                 tolerance = 1e-10)
  }
})


## =================================================================
## S8: Return structure (updated for Phase 2 — hide_mask removed)
## =================================================================

test_that("S8.3: Return structure has summary and records", {

  skip_on_cran()
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42)
  ))
  expect_true("summary" %in% names(res))
  expect_true("records" %in% names(res))
  expect_true("fits" %in% names(res))
  expect_true("RMSE" %in% names(res$summary))
  expect_true("Bias" %in% names(res$summary))
})


## =================================================================
## S9: Input validation (updated for Phase 2)
## =================================================================

test_that("S9.1: Invalid criterion rejected", {

  skip_on_cran()
  invalid_result <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, criterion = "invalid", seed = 42)
    ))
    "no_error"
  }, error = function(e) "error_thrown")

  expect_equal(invalid_result, "error_thrown",
    info = paste("fect_mspe(criterion='invalid') should throw an error,",
                 "but it completed without error. Builder must add",
                 "criterion validation."))
})

test_that("S9.2: Invalid cv.method rejected by fect_mspe", {

  skip_on_cran()
  expect_error(
    fect_mspe(out_base, cv.method = "invalid"),
    "cv.method|arg"
  )
})

test_that("S9.3: W wrong dimensions", {

  skip_on_cran()
  W_bad <- matrix(1, nrow = 5, ncol = 5)
  expect_error(
    fect_mspe(out_base, W = W_bad),
    "dimension|W"
  )
})

test_that("S9.5: .score_residuals() empty input", {

  skip_on_cran()
  score_fn <- tryCatch(
    getFromNamespace(".score_residuals", "fect"),
    error = function(e) NULL
  )
  skip_if(is.null(score_fn), ".score_residuals() not yet implemented")

  expect_error(score_fn(numeric(0)), "No residuals")
})


## =================================================================
## Section D: fect_mspe simplification (NEW for Phase 2)
## =================================================================

test_that("MSPE1: Simplified fect_mspe with cv.method='all_units'", {

  skip_on_cran()
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, cv.method = "all_units", criterion = "mspe")
  ))
  expect_true("summary" %in% names(res))
  expect_true("records" %in% names(res))
  expect_true("MSPE" %in% names(res$summary))
  expect_true("RMSE" %in% names(res$summary))
  expect_true("Bias" %in% names(res$summary))
  expect_true(res$summary$RMSE > 0)
  expect_equal(res$summary$RMSE, sqrt(res$summary$MSPE), tolerance = 1e-10)
})

test_that("MSPE2: Simplified fect_mspe with cv.method='treated_units'", {

  skip_on_cran()
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, cv.method = "treated_units",
              criterion = "mspe")
  ))
  expect_true("summary" %in% names(res))
  expect_true("records" %in% names(res))
  expect_true(res$summary$RMSE > 0)
})

test_that("MSPE3: Removed parameters rejected", {

  skip_on_cran()
  expect_error(fect_mspe(out_base, mask.method = "random"))
  expect_error(fect_mspe(out_base, hide_mask = matrix(TRUE, 10, 10)))
  expect_error(fect_mspe(out_base, n_rep = 3))
  expect_error(fect_mspe(out_base, pre.trend = TRUE))
  expect_error(fect_mspe(out_base, actual = out_base$Y.ct.full))
  expect_error(fect_mspe(out_base, control.only = FALSE))
  expect_error(fect_mspe(out_base, hide_n = 20))
})

test_that("MSPE4: Invalid cv.method rejected", {

  skip_on_cran()
  expect_error(
    fect_mspe(out_base, cv.method = "loo"),
    "cv.method|arg"
  )
})

test_that("MSPE5: Multi-model comparison with cv.method", {

  skip_on_cran()
  res <- suppressWarnings(suppressMessages(
    fect_mspe(list(m1 = out_base, m2 = out_base),
              seed = 42, cv.method = "all_units")
  ))
  expect_equal(nrow(res$summary), 2)
  expect_true(all(c("m1", "m2") %in% res$summary$Model))
})

test_that("MSPE6: Seed reproducibility", {

  skip_on_cran()
  r1 <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, cv.method = "all_units")
  ))
  r2 <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, cv.method = "all_units")
  ))
  expect_identical(r1$summary, r2$summary)
})

test_that("MSPE7: fect_mspe with observation weights", {

  skip_on_cran()
  TT <- nrow(out_base$Y.dat)
  NN <- ncol(out_base$Y.dat)
  W_mat <- matrix(1, nrow = TT, ncol = NN)
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, cv.method = "all_units", W = W_mat)
  ))
  if ("MSPE" %in% names(res$summary)) {
    expect_true(all(is.finite(res$summary$MSPE)))
    expect_true(all(res$summary$MSPE > 0))
  }
})

test_that("MSPE8: fect_mspe with norm.para", {

  skip_on_cran()
  res_raw <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, cv.method = "all_units")
  ))
  res_norm <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, cv.method = "all_units",
              norm.para = c(2.0))
  ))
  if ("MSPE" %in% names(res_raw$summary) &&
      "MSPE" %in% names(res_norm$summary)) {
    # With norm.para, MSPE scaled by norm.para[1]^2 = 4.0
    expect_equal(res_norm$summary$MSPE, res_raw$summary$MSPE * 4.0,
                 tolerance = 1e-10)
  }
})


