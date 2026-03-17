## ---------------------------------------------------------------
## Tests for score unification: .score_residuals(), fect_cv
## regression, fect_mspe criterion/masking/weights extensions.
##
## Follows test-spec.md for REQ-score-unify-001.
## ---------------------------------------------------------------

## Shared fixture — fitted once, reused across blocks.
suppressWarnings(data("simdata", package = "fect"))

out_base <- suppressWarnings(suppressMessages(
  fect::fect(
    Y ~ D + X1 + X2,
    data    = simdata,
    index   = c("id", "time"),
    method  = "ife",
    r       = 2,
    CV      = FALSE,
    se      = FALSE,
    parallel = FALSE
  )
))

## =================================================================
## S1: .score_residuals() unit tests with known inputs
## =================================================================

test_that("S1.1: Basic unweighted scoring", {
  # .score_residuals is an internal (non-exported) function.
  # It will be created by builder in R/score.R.
  # Skip gracefully if not yet available.
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
  score_fn <- tryCatch(
    getFromNamespace(".score_residuals", "fect"),
    error = function(e) NULL
  )
  skip_if(is.null(score_fn), ".score_residuals() not yet implemented")

  expect_error(score_fn(numeric(0)), "No residuals")
})


## =================================================================
## S2: fect_cv regression test (before vs after refactor)
## =================================================================

test_that("S2.1: IFE method CV - r.cv and CV.out snapshot", {
  # Run CV with IFE to check that the refactored code still selects

  # the same r.cv and produces comparable MSPE values.
  cv_out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      CV        = TRUE,
      r         = c(0, 3),
      criterion = "mspe",
      se        = FALSE,
      parallel  = FALSE
    )
  ))
  # r.cv should be a non-negative integer in [0, 3]
  expect_true(cv_out$r.cv >= 0 && cv_out$r.cv <= 3)

  # CV.out should exist and have MSPE column with finite positive values
  expect_true(!is.null(cv_out$CV.out))
  mspe_col <- cv_out$CV.out[, "MSPE"]
  # At least some entries should be less than 1e20 (the init value)
  expect_true(any(mspe_col < 1e19))
  expect_true(all(is.finite(mspe_col[mspe_col < 1e19])))
  expect_true(all(mspe_col[mspe_col < 1e19] > 0))
})

test_that("S2.2: MC method CV - lambda.cv selection", {
  cv_mc <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2,
      data      = simdata,
      index     = c("id", "time"),
      method    = "mc",
      CV        = TRUE,
      criterion = "mspe",
      se        = FALSE,
      parallel  = FALSE
    )
  ))
  # lambda.cv should be selected
  expect_true(!is.null(cv_mc$lambda.cv) || !is.null(cv_mc$r.cv))
  expect_true(!is.null(cv_mc$CV.out))
})

test_that("S2.3: GMoment column correctly populated (IFE)", {
  # After refactor, CV.out[,"GMoment"] should contain actual GMoment
  # values, not copies of MSPTATT.
  cv_out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      CV        = TRUE,
      r         = c(0, 3),
      criterion = "mspe",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  # GMoment values should be finite and positive where computed
  gm_col <- cv_out$CV.out[, "GMoment"]
  computed <- gm_col[gm_col < 1e19]
  if (length(computed) > 0) {
    expect_true(all(is.finite(computed)))
    expect_true(all(computed > 0))

    # GMoment should generally differ from MSPTATT
    msptatt_col <- cv_out$CV.out[, "MSPTATT"]
    msptatt_computed <- msptatt_col[gm_col < 1e19]
    if (length(msptatt_computed) > 0) {
      # At least some values should differ (bug fix verification)
      # Note: they could coincidentally match for some r, so we just
      # check that they are not ALL identical
      if (length(computed) > 1) {
        expect_false(
          all(abs(computed - msptatt_computed) < 1e-15),
          info = "GMoment should not be identical to MSPTATT for all r values"
        )
      }
    }
  }
})


## =================================================================
## S3: fect_mspe criterion support
## =================================================================

test_that("S3.1: Default criterion='mspe' matches old RMSE", {
  # Skip if criterion param not yet implemented
  has_criterion <- tryCatch({
    res <- suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 123, hide_n = 20, n_rep = 1,
                criterion = "mspe")
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_criterion, "criterion parameter not yet implemented in fect_mspe")

  res_old <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 123, hide_n = 20, n_rep = 1)
  ))
  res_new <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 123, hide_n = 20, n_rep = 1,
              criterion = "mspe")
  ))

  # Backward compat: RMSE should match
  expect_equal(res_new$summary$RMSE, res_old$summary$RMSE, tolerance = 1e-10)

  # RMSE = sqrt(MSPE) invariant (P1)
  if ("MSPE" %in% names(res_new$summary)) {
    expect_equal(res_new$summary$RMSE, sqrt(res_new$summary$MSPE),
                 tolerance = 1e-10)
  }
})

test_that("S3.2: All 7 criteria produce finite, positive scores", {
  has_criterion <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 123, hide_n = 20,
                criterion = "mspe")
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_criterion, "criterion parameter not yet implemented in fect_mspe")

  for (crit in c("mspe", "wmspe", "gmspe", "wgmspe", "mad",
                 "moment", "gmoment")) {
    res <- suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 123, hide_n = 20, criterion = crit)
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
## S4: fect_mspe with mask.method="cv.sample"
## =================================================================

test_that("S4.1: cv.sample masking runs without error", {
  has_mask <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 42, mask.method = "cv.sample",
                k = 3, cv.prop = 0.1, cv.treat = TRUE,
                cv.nobs = 3, cv.donut = 1, min.T0 = 3,
                criterion = "mspe")
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_mask, "mask.method parameter not yet implemented in fect_mspe")

  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, mask.method = "cv.sample",
              k = 3, cv.prop = 0.1, cv.treat = TRUE,
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

test_that("S4.2: cv.sample differs from random masking", {
  has_mask <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 42, mask.method = "cv.sample",
                k = 3, criterion = "mspe")
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_mask, "mask.method not yet implemented")

  res_random <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20,
              mask.method = "random", criterion = "mspe")
  ))
  res_cv <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, mask.method = "cv.sample",
              k = 3, criterion = "mspe")
  ))
  # Scores should differ since masking strategies are different
  expect_false(identical(res_random$summary$RMSE, res_cv$summary$RMSE))
})

test_that("S4.3: cv.sample with k=1 (random fold)", {
  has_mask <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 42, mask.method = "cv.sample",
                k = 1, criterion = "mspe")
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_mask, "mask.method not yet implemented")

  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, mask.method = "cv.sample",
              k = 1, criterion = "mspe")
  ))
  expect_true(res$summary$RMSE > 0)
})


## =================================================================
## S5: fect_mspe with actual and control.only (fect_mspe_sim replacement)
## =================================================================

test_that("S5.1: actual parameter with ground-truth matrix", {
  has_actual <- tryCatch({
    Y0_true <- out_base$Y.ct.full
    Y0_perturbed <- Y0_true + matrix(
      rnorm(length(Y0_true), sd = 0.1),
      nrow = nrow(Y0_true)
    )
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 42, hide_n = 20,
                actual = Y0_perturbed, control.only = TRUE,
                criterion = "mspe")
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_actual, "actual parameter not yet implemented")

  set.seed(99)
  Y0_true <- out_base$Y.ct.full
  Y0_perturbed <- Y0_true + matrix(
    rnorm(length(Y0_true), sd = 0.1),
    nrow = nrow(Y0_true)
  )
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20,
              actual = Y0_perturbed, control.only = TRUE,
              criterion = "mspe")
  ))
  if ("MSPE" %in% names(res$summary)) {
    expect_true(is.finite(res$summary$MSPE))
    expect_true(res$summary$MSPE > 0)
  }
})

test_that("S5.2: control.only=FALSE masks treated cells", {
  has_actual <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 42, hide_n = 20,
                actual = out_base$Y.ct.full,
                control.only = FALSE, criterion = "mspe")
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_actual, "actual/control.only not yet implemented")

  Y0_true <- out_base$Y.ct.full
  res_all <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20,
              actual = Y0_true, control.only = FALSE,
              criterion = "mspe")
  ))
  if ("MSPE" %in% names(res_all$summary)) {
    expect_true(is.finite(res_all$summary$MSPE))
  }
})

test_that("S5.3: When actual equals Y.ct.full, results are finite", {
  has_actual <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 42, hide_n = 20,
                actual = out_base$Y.ct.full,
                control.only = TRUE, criterion = "mspe")
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_actual, "actual parameter not yet implemented")

  Y0_exact <- out_base$Y.ct.full
  res_exact <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20,
              actual = Y0_exact, control.only = TRUE,
              criterion = "mspe")
  ))
  if ("MSPE" %in% names(res_exact$summary)) {
    expect_true(is.finite(res_exact$summary$MSPE))
  }
})


## =================================================================
## S6: fect_mspe with observation weights (W)
## =================================================================

test_that("S6.1: W parameter produces different scores than unweighted", {
  has_w <- tryCatch({
    TT <- nrow(out_base$Y.dat)
    NN <- ncol(out_base$Y.dat)
    W_mat <- matrix(runif(TT * NN, 0.5, 1.5), nrow = TT, ncol = NN)
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 42, hide_n = 20,
                criterion = "mspe", W = W_mat)
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_w, "W parameter not yet implemented in fect_mspe")

  TT <- nrow(out_base$Y.dat)
  NN <- ncol(out_base$Y.dat)
  set.seed(7)
  W_mat <- matrix(runif(TT * NN, 0.5, 1.5), nrow = TT, ncol = NN)
  res_unw <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20, criterion = "mspe",
              W = NULL)
  ))
  res_w <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20, criterion = "mspe",
              W = W_mat)
  ))
  expect_false(identical(res_unw$summary$RMSE, res_w$summary$RMSE))
})

test_that("S6.2: Uniform W equals unweighted", {
  has_w <- tryCatch({
    TT <- nrow(out_base$Y.dat)
    NN <- ncol(out_base$Y.dat)
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 42, hide_n = 20,
                criterion = "mspe",
                W = matrix(1, nrow = TT, ncol = NN))
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_w, "W parameter not yet implemented")

  TT <- nrow(out_base$Y.dat)
  NN <- ncol(out_base$Y.dat)
  W_uniform <- matrix(1, nrow = TT, ncol = NN)
  res_u <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20,
              criterion = "mspe", W = W_uniform)
  ))
  res_n <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20,
              criterion = "mspe", W = NULL)
  ))
  if ("MSPE" %in% names(res_u$summary) && "MSPE" %in% names(res_n$summary)) {
    expect_equal(res_u$summary$MSPE, res_n$summary$MSPE, tolerance = 1e-10)
  }
})


## =================================================================
## S7: fect_mspe with norm.para
## =================================================================

test_that("S7.1: norm.para scales scores", {
  has_np <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 42, hide_n = 20,
                criterion = "mspe", norm.para = c(2.0, 0.0))
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_np, "norm.para not yet implemented in fect_mspe")

  np <- c(2.0, 0.0)
  res_raw <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20, criterion = "mspe")
  ))
  res_norm <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20, criterion = "mspe",
              norm.para = np)
  ))
  if ("MSPE" %in% names(res_raw$summary) &&
      "MSPE" %in% names(res_norm$summary)) {
    expect_equal(res_norm$summary$MSPE, res_raw$summary$MSPE * 4.0,
                 tolerance = 1e-10)
  }
})


## =================================================================
## S8: Backward compatibility
## =================================================================

test_that("S8.1: Old-style pre.trend=TRUE still works", {
  has_mask <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 42, hide_n = 20,
                mask.method = "pre.trend", pre.trend.n = 2)
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_mask, "mask.method not yet implemented")

  res_old <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20,
              pre.trend = TRUE, pre.trend.n = 2)
  ))
  res_new <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20,
              mask.method = "pre.trend", pre.trend.n = 2)
  ))
  expect_equal(res_old$summary$RMSE, res_new$summary$RMSE, tolerance = 1e-10)
})

test_that("S8.2: Old-style pre.trend=FALSE still works", {
  has_mask <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 42, hide_n = 20,
                mask.method = "random")
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_mask, "mask.method not yet implemented")

  res_old <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20,
              pre.trend = FALSE)
  ))
  res_new <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20,
              mask.method = "random")
  ))
  expect_equal(res_old$summary$RMSE, res_new$summary$RMSE, tolerance = 1e-10)
})

test_that("S8.3: Old return structure preserved", {
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20)
  ))
  expect_true("summary" %in% names(res))
  expect_true("records" %in% names(res))
  expect_true("hide_mask" %in% names(res))
  expect_true("fits" %in% names(res))
  expect_true("RMSE" %in% names(res$summary))
  expect_true("Bias" %in% names(res$summary))
})


## =================================================================
## S9: Input validation
## =================================================================

test_that("S9.1: Invalid criterion rejected", {
  # Builder must add input validation: an invalid criterion value
  # should produce an error, not silently accept it.
  # This test will FAIL if fect_mspe silently accepts invalid criteria.
  invalid_result <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, criterion = "invalid", seed = 42, hide_n = 5)
    ))
    "no_error"
  }, error = function(e) "error_thrown")

  # If no error is thrown, that means criterion validation is missing.
  # The test spec requires that invalid criterion values produce an error.
  expect_equal(invalid_result, "error_thrown",
    info = paste("fect_mspe(criterion='invalid') should throw an error,",
                 "but it completed without error. Builder must add",
                 "criterion validation."))
})

test_that("S9.2: Invalid mask.method", {
  has_mask <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, mask.method = "random")
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_mask, "mask.method not yet implemented")

  expect_error(
    fect_mspe(out_base, mask.method = "invalid"),
    "mask.method|arg"
  )
})

test_that("S9.3: W wrong dimensions", {
  has_w <- tryCatch({
    TT <- nrow(out_base$Y.dat)
    NN <- ncol(out_base$Y.dat)
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, criterion = "mspe",
                W = matrix(1, nrow = TT, ncol = NN))
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_w, "W parameter not yet implemented")

  W_bad <- matrix(1, nrow = 5, ncol = 5)
  expect_error(
    fect_mspe(out_base, W = W_bad),
    "dimension|W"
  )
})

test_that("S9.4: actual wrong dimensions", {
  has_actual <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, actual = out_base$Y.ct.full,
                control.only = TRUE, criterion = "mspe")
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_actual, "actual parameter not yet implemented")

  actual_bad <- matrix(1, nrow = 5, ncol = 5)
  expect_error(
    fect_mspe(out_base, actual = actual_bad),
    "dimension|actual"
  )
})

test_that("S9.5: .score_residuals() empty input", {
  score_fn <- tryCatch(
    getFromNamespace(".score_residuals", "fect"),
    error = function(e) NULL
  )
  skip_if(is.null(score_fn), ".score_residuals() not yet implemented")

  expect_error(score_fn(numeric(0)), "No residuals")
})


## =================================================================
## Property-based invariants
## =================================================================

test_that("P1: RMSE = sqrt(MSPE) invariant", {
  has_criterion <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 42, hide_n = 20, criterion = "mspe")
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_criterion, "criterion not yet implemented")

  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20, criterion = "mspe")
  ))
  if ("MSPE" %in% names(res$summary)) {
    expect_equal(res$summary$RMSE, sqrt(res$summary$MSPE), tolerance = 1e-10)
  }
})

test_that("P2: MSPE >= 0", {
  has_criterion <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 42, hide_n = 20, criterion = "mspe")
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_criterion, "criterion not yet implemented")

  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20, criterion = "mspe")
  ))
  if ("MSPE" %in% names(res$summary)) {
    expect_true(all(res$summary$MSPE >= 0))
  }
})

test_that("P5: Moment = 0 when all per-group mean residuals are 0", {
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
  r1 <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20)
  ))
  r2 <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 20)
  ))
  expect_identical(r1$summary, r2$summary)
})


## =================================================================
## Edge cases
## =================================================================

test_that("E3: Single model in fect_mspe (not a list)", {
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 10)
  ))
  expect_equal(nrow(res$summary), 1)
  expect_true(res$summary$RMSE > 0)
})

test_that("E4: Multiple models in fect_mspe", {
  multi <- suppressWarnings(suppressMessages(
    fect_mspe(list(m1 = out_base, m2 = out_base),
              seed = 42, hide_n = 10)
  ))
  expect_equal(nrow(multi$summary), 2)
  expect_true(all(c("m1", "m2") %in% multi$summary$Model))
})

test_that("E5: n_rep > 1 with standard masking", {
  has_criterion <- tryCatch({
    suppressWarnings(suppressMessages(
      fect_mspe(out_base, seed = 42, hide_n = 10,
                n_rep = 3, criterion = "mspe")
    ))
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_criterion, "criterion not yet implemented")

  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, hide_n = 10,
              n_rep = 3, criterion = "mspe")
  ))
  # 3 reps x 1 model = 3 records
  expect_equal(nrow(res$records), 3)
  # Summary should aggregate across reps
  expect_equal(nrow(res$summary), 1)
})
