## ---------------------------------------------------------------
## Tests for .score_residuals() unit behavior + property-based
## invariants + edge cases.
##
## Originally Section S1 + Property + Edge of test-score-unify.R;
## split out 2026-05-03 for progress visibility under reporter = "summary".
## Shared fixtures (make_factor_data, ntdata, out_base) live in
## helper-score-unify.R.
## ---------------------------------------------------------------

## =================================================================
## S1: .score_residuals() unit tests with known inputs
## =================================================================

test_that("S1.1: Basic unweighted scoring", {

  skip_on_cran()
  score_fn <- tryCatch(
    getFromNamespace(".score_residuals", "fect"),
    error = function(e) NULL
  )
  skip_if(is.null(score_fn), ".score_residuals() not yet implemented")

  resid <- c(1.0, -2.0, 3.0, -1.0, 2.0)
  result <- score_fn(resid)

  # MSPE = mean(resid^2) = (1+4+9+1+4)/5 = 3.8
  expect_equal(result[["MSPE"]], 3.8, tolerance = 1e-10)

  # GMSPE = exp(mean(log(resid^2)))
  expected_gmspe <- exp(mean(log(c(1, 4, 9, 1, 4))))
  expect_equal(result[["GMSPE"]], expected_gmspe, tolerance = 1e-10)

  # MAD = median(|e2 - median(e2)|) where e2=c(1,4,9,1,4), median=4
  # deviations = c(3,0,5,3,0), MAD = 3
  expect_equal(result[["MAD"]], 3.0, tolerance = 1e-10)

  # RMSE = sqrt(3.8)
  expect_equal(result[["RMSE"]], sqrt(3.8), tolerance = 1e-10)

  # Bias = mean(resid) = 0.6
  expect_equal(result[["Bias"]], 0.6, tolerance = 1e-10)

  # Moment and GMoment = NA (no time_index)
  expect_true(is.na(result[["Moment"]]))
  expect_true(is.na(result[["GMoment"]]))

  # WMSPE = MSPE with uniform weights = 3.8
  expect_equal(result[["WMSPE"]], 3.8, tolerance = 1e-10)

  # WGMSPE = GMSPE with uniform weights
  expect_equal(result[["WGMSPE"]], expected_gmspe, tolerance = 1e-10)
})

test_that("S1.2: With observation weights", {

  skip_on_cran()
  score_fn <- tryCatch(
    getFromNamespace(".score_residuals", "fect"),
    error = function(e) NULL
  )
  skip_if(is.null(score_fn), ".score_residuals() not yet implemented")

  resid <- c(1.0, -2.0, 3.0)
  obs_weights <- c(0.5, 1.0, 0.5)
  result <- score_fn(resid, obs_weights = obs_weights)

  # MSPE = (0.5*1 + 1.0*4 + 0.5*9) / (0.5+1.0+0.5) = 9/2 = 4.5
  expect_equal(result[["MSPE"]], 4.5, tolerance = 1e-10)

  # RMSE = sqrt(4.5)
  expect_equal(result[["RMSE"]], sqrt(4.5), tolerance = 1e-10)

  # Bias = mean(resid) = 2/3
  expect_equal(result[["Bias"]], 2 / 3, tolerance = 1e-10)
})

test_that("S1.3: With time_index and count_weights (Moment/GMoment)", {

  skip_on_cran()
  score_fn <- tryCatch(
    getFromNamespace(".score_residuals", "fect"),
    error = function(e) NULL
  )
  skip_if(is.null(score_fn), ".score_residuals() not yet implemented")

  resid <- c(1.0, -1.0, 2.0, -3.0, 0.5, -0.5)
  time_index <- c("-2", "-2", "-1", "-1", "Control", "Control")
  count_weights <- c("-2" = 1.5, "-1" = 2.0, "Control" = 1.0)

  result <- score_fn(resid, time_index = time_index,
                     count_weights = count_weights)

  # Moment:
  # resid_mean per group: "-2"->0, "-1"->-0.5, "Control"->0
  # abs: c(0, 0.5, 0)
  # Moment = (1.5*0 + 2.0*0.5 + 1.0*0) / (1.5+2.0+1.0) = 1.0/4.5
  expected_moment <- 1.0 / 4.5
  expect_equal(result[["Moment"]], expected_moment, tolerance = 1e-10)

  # GMoment:
  # geometric mean of abs(resid) per group:
  # "-2" -> exp((log(1)+log(1))/2) = 1.0
  # "-1" -> exp((log(2)+log(3))/2) = sqrt(6)
  # "Control" -> exp((log(0.5)+log(0.5))/2) = 0.5
  # GMoment = (1.5*1.0 + 2.0*sqrt(6) + 1.0*0.5) / (1.5+2.0+1.0)
  expected_gmoment <- (1.5 * 1.0 + 2.0 * sqrt(6) + 1.0 * 0.5) / 4.5
  expect_equal(result[["GMoment"]], expected_gmoment, tolerance = 1e-10)
})

test_that("S1.4: With norm.para", {

  skip_on_cran()
  score_fn <- tryCatch(
    getFromNamespace(".score_residuals", "fect"),
    error = function(e) NULL
  )
  skip_if(is.null(score_fn), ".score_residuals() not yet implemented")

  resid <- c(1.0, -2.0, 3.0, -1.0, 2.0)
  result <- score_fn(resid, norm.para = c(2.0))

  # All 7 scores multiplied by 4.0 (norm.para[1]^2)
  expect_equal(result[["MSPE"]], 3.8 * 4.0, tolerance = 1e-10)
  expect_equal(result[["RMSE"]], sqrt(3.8 * 4.0), tolerance = 1e-10)

  # Bias unchanged
  expect_equal(result[["Bias"]], 0.6, tolerance = 1e-10)
})

test_that("S1.5: Edge case - single residual", {

  skip_on_cran()
  score_fn <- tryCatch(
    getFromNamespace(".score_residuals", "fect"),
    error = function(e) NULL
  )
  skip_if(is.null(score_fn), ".score_residuals() not yet implemented")

  result <- score_fn(c(2.5))

  expect_equal(result[["MSPE"]], 6.25, tolerance = 1e-10)
  expect_equal(result[["RMSE"]], 2.5, tolerance = 1e-10)
  expect_equal(result[["Bias"]], 2.5, tolerance = 1e-10)
  expect_equal(result[["MAD"]], 0.0, tolerance = 1e-10)
  expect_equal(result[["GMSPE"]], 6.25, tolerance = 1e-10)
})

test_that("S1.6: Edge case - zero residual", {

  skip_on_cran()
  score_fn <- tryCatch(
    getFromNamespace(".score_residuals", "fect"),
    error = function(e) NULL
  )
  skip_if(is.null(score_fn), ".score_residuals() not yet implemented")

  resid <- c(0.0, 1.0, -1.0)
  result <- score_fn(resid)

  # MSPE = 2/3
  expect_equal(result[["MSPE"]], 2 / 3, tolerance = 1e-10)

  # GMSPE = exp(mean(log(c(0,1,1)))) = exp(-Inf) = 0
  expect_true(result[["GMSPE"]] <= 1e-300 || !is.finite(log(result[["GMSPE"]])))

  # WGMSPE: zero filtered out, only c(1,1) remain. WGMSPE = 1.0
  expect_equal(result[["WGMSPE"]], 1.0, tolerance = 1e-10)
})

test_that("S1.7: Edge case - empty residuals", {

  skip_on_cran()
  score_fn <- tryCatch(
    getFromNamespace(".score_residuals", "fect"),
    error = function(e) NULL
  )
  skip_if(is.null(score_fn), ".score_residuals() not yet implemented")

  expect_error(score_fn(numeric(0)), "No residuals")
})


## =================================================================
## Property-based invariants (updated for Phase 2)
## =================================================================

test_that("P1: RMSE = sqrt(MSPE) invariant", {

  skip_on_cran()
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, criterion = "mspe")
  ))
  if ("MSPE" %in% names(res$summary)) {
    expect_equal(res$summary$RMSE, sqrt(res$summary$MSPE), tolerance = 1e-10)
  }
})

test_that("P2: MSPE >= 0", {

  skip_on_cran()
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, criterion = "mspe")
  ))
  if ("MSPE" %in% names(res$summary)) {
    expect_true(all(res$summary$MSPE >= 0))
  }
})

test_that("P5: Moment = 0 when all per-group mean residuals are 0", {

  skip_on_cran()
  score_fn <- tryCatch(
    getFromNamespace(".score_residuals", "fect"),
    error = function(e) NULL
  )
  skip_if(is.null(score_fn), ".score_residuals() not yet implemented")

  # Construct residuals where each group mean is exactly 0
  resid <- c(1, -1, 2, -2)
  time_index <- c("A", "A", "B", "B")
  count_weights <- c("A" = 1.0, "B" = 1.0)

  result <- score_fn(resid, time_index = time_index,
                     count_weights = count_weights)
  expect_equal(result[["Moment"]], 0.0, tolerance = 1e-10)
})

test_that("P7: Seed reproducibility", {

  skip_on_cran()
  r1 <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42)
  ))
  r2 <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42)
  ))
  expect_identical(r1$summary, r2$summary)
})


## =================================================================
## Edge cases (updated for Phase 2)
## =================================================================

test_that("E3: Single model in fect_mspe (not a list)", {

  skip_on_cran()
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42)
  ))
  expect_equal(nrow(res$summary), 1)
  expect_true(res$summary$RMSE > 0)
})

test_that("E4: Multiple models in fect_mspe", {

  skip_on_cran()
  multi <- suppressWarnings(suppressMessages(
    fect_mspe(list(m1 = out_base, m2 = out_base),
              seed = 42)
  ))
  expect_equal(nrow(multi$summary), 2)
  expect_true(all(c("m1", "m2") %in% multi$summary$Model))
})


