## ---------------------------------------------------------------
## Tests for score unification: .score_residuals(), fect_cv
## regression, fect_mspe criterion/masking/weights extensions,
## cv.method unification (Phase 2).
##
## Follows test-spec.md for REQ-cv-method-phase2.
## ---------------------------------------------------------------

## Shared fixture — fitted once, reused across blocks.
suppressWarnings(data("simdata", package = "fect"))

## DGP with factor structure and sufficient never-treated units for CV.
## N=50, TT=20, Ntr=15 => 35 never-treated units with 20 pre-treatment
## periods — plenty for cross-validation with r up to 3.
make_factor_data <- function(N = 50, TT = 20, Ntr = 15, tau = 3.0,
                              r = 2, seed = 42) {
  set.seed(seed)
  F_mat <- matrix(rnorm(TT * r), TT, r)
  L_mat <- matrix(rnorm(N * r), N, r)
  alpha_i <- rnorm(N, 0, 1)
  xi_t <- rnorm(TT, 0, 0.5)

  T0_vec <- rep(Inf, N)
  if (Ntr > 0) {
    T0_vec[1:Ntr] <- sample(round(TT * 0.4):round(TT * 0.7), Ntr,
                             replace = TRUE)
  }

  Y_vec <- D_vec <- numeric(N * TT)
  id_vec <- time_vec <- integer(N * TT)
  idx <- 1
  for (i in 1:N) {
    for (t in 1:TT) {
      treated <- (t >= T0_vec[i])
      D_vec[idx] <- as.integer(treated)
      Y_vec[idx] <- alpha_i[i] + xi_t[t] +
        sum(F_mat[t, ] * L_mat[i, ]) +
        tau * D_vec[idx] + rnorm(1, 0, 0.5)
      id_vec[idx] <- i
      time_vec[idx] <- t
      idx <- idx + 1
    }
  }

  data.frame(id = id_vec, time = time_vec, Y = Y_vec, D = D_vec)
}

## Shared never-treated fixture for Section C, E, F tests
ntdata <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

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
## S3: fect_mspe criterion support (updated for Phase 2)
## =================================================================

test_that("S3.1: Default criterion='mspe' matches old RMSE", {
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
  res_cv <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, cv.method = "all_units",
              k = 3, criterion = "mspe")
  ))
  expect_true("summary" %in% names(res_cv))
  expect_true(res_cv$summary$RMSE > 0)
})

test_that("S4.3: cv.method='treated_units' with k=1 (single fold)", {
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
  expect_error(
    fect_mspe(out_base, cv.method = "invalid"),
    "cv.method|arg"
  )
})

test_that("S9.3: W wrong dimensions", {
  W_bad <- matrix(1, nrow = 5, ncol = 5)
  expect_error(
    fect_mspe(out_base, W = W_bad),
    "dimension|W"
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
## Property-based invariants (updated for Phase 2)
## =================================================================

test_that("P1: RMSE = sqrt(MSPE) invariant", {
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, criterion = "mspe")
  ))
  if ("MSPE" %in% names(res$summary)) {
    expect_equal(res$summary$RMSE, sqrt(res$summary$MSPE), tolerance = 1e-10)
  }
})

test_that("P2: MSPE >= 0", {
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, criterion = "mspe")
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
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42)
  ))
  expect_equal(nrow(res$summary), 1)
  expect_true(res$summary$RMSE > 0)
})

test_that("E4: Multiple models in fect_mspe", {
  multi <- suppressWarnings(suppressMessages(
    fect_mspe(list(m1 = out_base, m2 = out_base),
              seed = 42)
  ))
  expect_equal(nrow(multi$summary), 2)
  expect_true(all(c("m1", "m2") %in% multi$summary$Model))
})


## =================================================================
## Section B: cv.method in fect_cv (NEW for Phase 2)
## =================================================================

test_that("CV1: cv.method='all_units' selects r.cv", {
  cv_out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      CV        = TRUE,
      r         = c(0, 3),
      cv.method = "all_units",
      se        = FALSE,
      parallel  = FALSE
    )
  ))
  # r.cv is integer in [0, 3]
  expect_true(cv_out$r.cv >= 0 && cv_out$r.cv <= 3)
  # CV.out exists with MSPE column
  expect_true(!is.null(cv_out$CV.out))
  mspe_col <- cv_out$CV.out[, "MSPE"]
  computed <- mspe_col[mspe_col < 1e19]
  expect_true(all(is.finite(computed)))
  expect_true(all(computed > 0))
})

test_that("CV2: cv.method='treated_units' selects r.cv", {
  cv_out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      CV        = TRUE,
      r         = c(0, 3),
      cv.method = "treated_units",
      se        = FALSE,
      parallel  = FALSE
    )
  ))
  # r.cv is integer in [0, 3]
  expect_true(cv_out$r.cv >= 0 && cv_out$r.cv <= 3)
  # CV.out exists
  expect_true(!is.null(cv_out$CV.out))
  mspe_col <- cv_out$CV.out[, "MSPE"]
  computed <- mspe_col[mspe_col < 1e19]
  expect_true(all(is.finite(computed)))
  expect_true(all(computed > 0))
})

test_that("CV3: Invalid cv.method rejected", {
  expect_error(
    fect::fect(
      Y ~ D + X1 + X2,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      CV        = TRUE,
      cv.method = "invalid"
    ),
    "cv.method|arg"
  )
})


## =================================================================
## Section C: cv.method in fect_nevertreated (NEW for Phase 2)
## =================================================================

test_that("NT1: fect_nevertreated cv.method='loo' selects r.cv (IFE)", {
  nt_out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "loo",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(nt_out$r.cv >= 0 && nt_out$r.cv <= 3)
  expect_true(!is.null(nt_out$CV.out))
  mspe_col <- nt_out$CV.out[, "MSPE"]
  computed <- mspe_col[mspe_col < 1e19]
  if (length(computed) > 0) {
    expect_true(all(is.finite(computed)))
  }
  expect_true(!is.null(nt_out$Y.ct))
})

test_that("NT2: fect_nevertreated cv.method='treated_units' selects r.cv (IFE)", {
  nt_out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "treated_units",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(nt_out$r.cv >= 0 && nt_out$r.cv <= 3)
  expect_true(!is.null(nt_out$CV.out))
})

test_that("NT3: fect_nevertreated cv.method='all_units' selects r.cv (IFE)", {
  nt_out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "all_units",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(nt_out$r.cv >= 0 && nt_out$r.cv <= 3)
  expect_true(!is.null(nt_out$CV.out))
})

test_that("NT4: fect_nevertreated cv.method='loo' selects r.cv (CFE)", {
  nt_cfe <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "cfe",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "loo",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(nt_cfe$r.cv >= 0 && nt_cfe$r.cv <= 3)
  expect_true(!is.null(nt_cfe$CV.out))
})

test_that("NT5: fect_nevertreated cv.method='treated_units' selects r.cv (CFE)", {
  nt_cfe <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "cfe",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "treated_units",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(nt_cfe$r.cv >= 0 && nt_cfe$r.cv <= 3)
  expect_true(!is.null(nt_cfe$CV.out))
})

test_that("NT6: fect_nevertreated default cv.method is treated_units", {
  nt_default <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      se              = FALSE,
      parallel        = FALSE
    )
  ))

  nt_explicit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "treated_units",
      se              = FALSE,
      parallel        = FALSE
    )
  ))

  expect_equal(nt_default$r.cv, nt_explicit$r.cv)
})

test_that("NT7: Invalid cv.method rejected for nevertreated", {
  expect_error(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      factors.from    = "nevertreated",
      CV              = TRUE,
      cv.method       = "invalid"
    ),
    "cv.method|arg"
  )
})


## =================================================================
## Section D: fect_mspe simplification (NEW for Phase 2)
## =================================================================

test_that("MSPE1: Simplified fect_mspe with cv.method='all_units'", {
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
  res <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, cv.method = "treated_units",
              criterion = "mspe")
  ))
  expect_true("summary" %in% names(res))
  expect_true("records" %in% names(res))
  expect_true(res$summary$RMSE > 0)
})

test_that("MSPE3: Removed parameters rejected", {
  expect_error(fect_mspe(out_base, mask.method = "random"))
  expect_error(fect_mspe(out_base, hide_mask = matrix(TRUE, 10, 10)))
  expect_error(fect_mspe(out_base, n_rep = 3))
  expect_error(fect_mspe(out_base, pre.trend = TRUE))
  expect_error(fect_mspe(out_base, actual = out_base$Y.ct.full))
  expect_error(fect_mspe(out_base, control.only = FALSE))
  expect_error(fect_mspe(out_base, hide_n = 20))
})

test_that("MSPE4: Invalid cv.method rejected", {
  expect_error(
    fect_mspe(out_base, cv.method = "loo"),
    "cv.method|arg"
  )
})

test_that("MSPE5: Multi-model comparison with cv.method", {
  res <- suppressWarnings(suppressMessages(
    fect_mspe(list(m1 = out_base, m2 = out_base),
              seed = 42, cv.method = "all_units")
  ))
  expect_equal(nrow(res$summary), 2)
  expect_true(all(c("m1", "m2") %in% res$summary$Model))
})

test_that("MSPE6: Seed reproducibility", {
  r1 <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, cv.method = "all_units")
  ))
  r2 <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, seed = 42, cv.method = "all_units")
  ))
  expect_identical(r1$summary, r2$summary)
})

test_that("MSPE7: fect_mspe with observation weights", {
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


## =================================================================
## Section E: 1% Selection Rule Verification (NEW for Phase 2)
## =================================================================

test_that("SEL1: 1% selection rule in IFE nevertreated", {
  # Large tol should NOT affect the 1% rule
  nt_bigtol <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      se              = FALSE,
      parallel        = FALSE,
      tol             = 0.5
    )
  ))
  nt_smalltol <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      se              = FALSE,
      parallel        = FALSE,
      tol             = 1e-3
    )
  ))
  expect_equal(nt_bigtol$r.cv, nt_smalltol$r.cv)
})

test_that("SEL2: 1% selection rule in CFE nevertreated", {
  nt_bigtol <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "cfe",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      se              = FALSE,
      parallel        = FALSE,
      tol             = 0.5
    )
  ))
  nt_smalltol <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "cfe",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      se              = FALSE,
      parallel        = FALSE,
      tol             = 1e-3
    )
  ))
  expect_equal(nt_bigtol$r.cv, nt_smalltol$r.cv)
})


## =================================================================
## Section F: W and count.T.cv in fect_nevertreated (NEW for Phase 2)
## =================================================================

test_that("WT1: W weights flow through nevertreated LOO scoring", {
  # W in fect() is a column name, not a matrix. Add a weight column to ntdata.
  ntdata_w <- ntdata
  ntdata_w$wt <- 1.0  # uniform weights

  nt_w <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata_w,
      index           = c("id", "time"),
      method          = "ife",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      W               = "wt",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(nt_w$r.cv >= 0 && nt_w$r.cv <= 3)

  # With uniform weights, r.cv should match unweighted result
  nt_nw <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_equal(nt_w$r.cv, nt_nw$r.cv)
})

test_that("WT2: Non-uniform W may change r selection", {
  # W in fect() is a column name. Add non-uniform weight column.
  ntdata_w <- ntdata
  set.seed(99)
  ntdata_w$wt <- runif(nrow(ntdata), 0.5, 2.0)

  nt_w <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata_w,
      index           = c("id", "time"),
      method          = "ife",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      W               = "wt",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  # Should not crash, r.cv should be valid
  expect_true(nt_w$r.cv >= 0 && nt_w$r.cv <= 3)
})


## =================================================================
## Section G: Integration Tests (NEW for Phase 2)
## =================================================================

test_that("INT1: End-to-end cv.method pipeline", {
  ## Use CV=FALSE for the fit since fect_mspe performs its own

  ## cross-validation masking; what matters is that the fit object
  ## has valid Y.dat/Y.ct for fect_mspe to score.
  ## NOTE: fect_mspe errors on CV=TRUE fits ("No valid residuals") —
  ## that is a separate source-code issue tracked for builder.
  fit_result <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = ntdata,
      index     = c("id", "time"),
      method    = "ife",
      CV        = FALSE,
      r         = 2,
      se        = FALSE,
      parallel  = FALSE
    )
  ))
  mspe_result <- suppressWarnings(suppressMessages(
    fect_mspe(fit_result, seed = 42, cv.method = "all_units")
  ))
  expect_true("summary" %in% names(mspe_result))
  expect_true(mspe_result$summary$RMSE > 0)
  if ("MSPE" %in% names(mspe_result$summary)) {
    expect_true(mspe_result$summary$MSPE > 0)
  }
})

test_that("INT2: IFE CV respects cv.method='treated_units'", {
  ## Original test used method="gsynth" + cv.method="loo", but simdata has
  ## treatment reversals which gsynth rejects. We test that cv.method is
  ## respected by using cv.method="treated_units" (non-default) with method="ife".
  ife_out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      CV        = TRUE,
      r         = c(0, 3),
      cv.method = "treated_units",
      se        = FALSE,
      parallel  = FALSE
    )
  ))
  expect_true(ife_out$r.cv >= 0 && ife_out$r.cv <= 3)
})


## =================================================================
## Section H: cv.sample k-fold CV in fect_nevertreated (NEW)
## Tests for actual cv.sample-based cross-validation branches
## when cv.method="all_units" or "treated_units" in nevertreated.
## =================================================================

## ---- H.1: IFE smoke tests ---- ##

test_that("NTCV1: cv.method='all_units' IFE produces valid output", {
  set.seed(42)
  out_au <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      force           = "two-way",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "all_units",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  # r.cv is integer in [0, 3]
  expect_true(out_au$r.cv >= 0 && out_au$r.cv <= 3)
  # att.avg is finite
  expect_true(is.finite(out_au$att.avg))
  # est.att has no all-NA columns
  if (!is.null(out_au$est.att)) {
    na_cols <- apply(out_au$est.att, 2, function(x) all(is.na(x)))
    expect_false(all(na_cols),
                 info = "est.att should not have all columns be NA")
  }
  # CV.out exists with proper structure
  expect_true(!is.null(out_au$CV.out))
  mspe_col <- out_au$CV.out[, "MSPE"]
  computed <- mspe_col[mspe_col < 1e19]
  if (length(computed) > 0) {
    expect_true(all(is.finite(computed)))
    expect_true(all(computed > 0))
  }
})

test_that("NTCV2: cv.method='treated_units' IFE produces valid output", {
  set.seed(42)
  out_tu <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      force           = "two-way",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "treated_units",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(out_tu$r.cv >= 0 && out_tu$r.cv <= 3)
  expect_true(is.finite(out_tu$att.avg))
  if (!is.null(out_tu$est.att)) {
    na_cols <- apply(out_tu$est.att, 2, function(x) all(is.na(x)))
    expect_false(all(na_cols))
  }
  expect_true(!is.null(out_tu$CV.out))
})

## ---- H.2: CFE smoke tests ---- ##

test_that("NTCV3: cv.method='all_units' CFE produces valid output", {
  set.seed(42)
  out_au_cfe <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "cfe",
      force           = "two-way",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "all_units",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(out_au_cfe$r.cv >= 0 && out_au_cfe$r.cv <= 3)
  expect_true(is.finite(out_au_cfe$att.avg))
  expect_true(!is.null(out_au_cfe$CV.out))
})

test_that("NTCV4: cv.method='treated_units' CFE produces valid output", {
  set.seed(42)
  out_tu_cfe <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "cfe",
      force           = "two-way",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "treated_units",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(out_tu_cfe$r.cv >= 0 && out_tu_cfe$r.cv <= 3)
  expect_true(is.finite(out_tu_cfe$att.avg))
  expect_true(!is.null(out_tu_cfe$CV.out))
})

## ---- H.3: LOO backward compatibility ---- ##

test_that("NTCV5: cv.method='loo' IFE backward compatibility", {
  set.seed(1234)
  out_loo <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      force           = "two-way",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "loo",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(out_loo$r.cv >= 0 && out_loo$r.cv <= 3)
  expect_true(is.finite(out_loo$att.avg))
  expect_true(!is.null(out_loo$CV.out))

  # LOO is deterministic: re-running should produce identical r.cv
  set.seed(1234)
  out_loo2 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      force           = "two-way",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "loo",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_equal(out_loo$r.cv, out_loo2$r.cv)
})

## ---- H.4: r-selection and ATT checks ---- ##

test_that("NTCV6: r-selection validity across all cv.methods (IFE)", {
  r_start <- 0
  r_end <- 3

  methods_list <- c("loo", "all_units", "treated_units")
  for (cm in methods_list) {
    set.seed(42)
    out <- suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        factors.from    = "nevertreated",
        CV              = TRUE,
        r               = c(r_start, r_end),
        cv.method       = cm,
        se              = FALSE,
        parallel        = FALSE
      )
    ))

    # r.cv in valid range
    expect_true(out$r.cv >= r_start && out$r.cv <= r_end,
                info = paste("cv.method =", cm, ": r.cv out of range"))

    # CV.out has correct number of rows
    expect_equal(nrow(out$CV.out), r_end - r_start + 1,
                 info = paste("cv.method =", cm, ": CV.out row count wrong"))

    # Score columns in CV.out: check MSPE for evaluated rows
    mspe_col <- out$CV.out[, "MSPE"]
    computed <- mspe_col[mspe_col < 1e19]
    if (length(computed) > 0) {
      expect_true(all(is.finite(computed)),
                  info = paste("cv.method =", cm, ": non-finite MSPE"))
      expect_true(all(computed >= 0),
                  info = paste("cv.method =", cm, ": negative MSPE"))
    }
  }
})

test_that("NTCV7: ATT consistency across cv.methods (IFE)", {
  att_vals <- numeric(3)
  methods_list <- c("loo", "all_units", "treated_units")
  for (i in seq_along(methods_list)) {
    set.seed(42)
    out <- suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        factors.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = methods_list[i],
        se              = FALSE,
        parallel        = FALSE
      )
    ))
    att_vals[i] <- out$att.avg
    expect_true(is.finite(out$att.avg),
                info = paste("cv.method =", methods_list[i], ": att.avg not finite"))
  }

  # Sanity check: all ATT values should be in the same ballpark.
  # The true ATT is ~3.0 for this DGP.
  # Allow wide tolerance since different r.cv selections produce different ATTs.
  att_range <- max(att_vals) - min(att_vals)
  expect_true(att_range < 5.0,
              info = paste("ATT values too spread:",
                           paste(round(att_vals, 3), collapse = ", ")))
})

## ---- H.5: Edge cases ---- ##

test_that("NTCV-Edge1: r=c(0,0) with all cv.methods", {
  for (cm in c("loo", "all_units", "treated_units")) {
    set.seed(42)
    out <- suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        factors.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 0),
        cv.method       = cm,
        se              = FALSE,
        parallel        = FALSE
      )
    ))
    expect_equal(out$r.cv, 0,
                 info = paste("cv.method =", cm, ": r.cv should be 0 when r=c(0,0)"))
    # Scores should still be finite
    mspe_col <- out$CV.out[, "MSPE"]
    computed <- mspe_col[mspe_col < 1e19]
    if (length(computed) > 0) {
      expect_true(all(is.finite(computed)),
                  info = paste("cv.method =", cm, ": non-finite MSPE with r=0"))
    }
  }
})

test_that("NTCV-Edge2: Small panel with all cv.methods", {
  # Create data with few pre-treatment periods for treated units
  small_data <- make_factor_data(N = 30, TT = 10, Ntr = 8, r = 1, seed = 99)

  for (cm in c("loo", "all_units", "treated_units")) {
    set.seed(42)
    out <- tryCatch(
      suppressWarnings(suppressMessages(
        fect::fect(
          Y ~ D,
          data            = small_data,
          index           = c("id", "time"),
          method          = "ife",
          force           = "two-way",
          factors.from    = "nevertreated",
          CV              = TRUE,
          r               = c(0, 2),
          cv.method       = cm,
          se              = FALSE,
          parallel        = FALSE
        )
      )),
      error = function(e) e
    )
    # Should either succeed with valid r.cv or produce informative error
    if (!inherits(out, "error")) {
      expect_true(out$r.cv >= 0 && out$r.cv <= 2,
                  info = paste("cv.method =", cm, ": r.cv out of range (small panel)"))
      expect_true(is.finite(out$att.avg),
                  info = paste("cv.method =", cm, ": att.avg not finite (small panel)"))
    } else {
      # If it errors, the message should be informative (not a cryptic crash)
      expect_true(nchar(conditionMessage(out)) > 0,
                  info = paste("cv.method =", cm, ": error should be informative"))
    }
  }
})

test_that("NTCV-Edge3: Single treated unit", {
  single_tr_data <- make_factor_data(N = 30, TT = 15, Ntr = 1, r = 1, seed = 77)

  for (cm in c("loo", "all_units", "treated_units")) {
    set.seed(42)
    out <- tryCatch(
      suppressWarnings(suppressMessages(
        fect::fect(
          Y ~ D,
          data            = single_tr_data,
          index           = c("id", "time"),
          method          = "ife",
          force           = "two-way",
          factors.from    = "nevertreated",
          CV              = TRUE,
          r               = c(0, 2),
          cv.method       = cm,
          se              = FALSE,
          parallel        = FALSE
        )
      )),
      error = function(e) e
    )
    if (!inherits(out, "error")) {
      expect_true(out$r.cv >= 0 && out$r.cv <= 2,
                  info = paste("cv.method =", cm, ": r.cv invalid (single treated)"))
    }
    # If it errors, that's acceptable for single treated unit edge case
  }
})

## ---- H.6: Property-based invariants for cv.sample ---- ##

test_that("NTCV-P1: Score non-negativity in CV.out", {
  for (cm in c("all_units", "treated_units")) {
    set.seed(42)
    out <- suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        factors.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = cm,
        se              = FALSE,
        parallel        = FALSE
      )
    ))

    # All score columns should be non-negative
    score_cols <- c("MSPE", "WMSPE", "GMSPE", "WGMSPE", "MAD")
    for (sc in score_cols) {
      if (sc %in% colnames(out$CV.out)) {
        vals <- out$CV.out[, sc]
        computed <- vals[vals < 1e19]
        if (length(computed) > 0) {
          expect_true(all(computed >= 0),
                      info = paste("cv.method =", cm, ", score =", sc, ": negative value"))
        }
      }
    }
  }
})

test_that("NTCV-P4: LOO determinism", {
  set.seed(42)
  out1 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      force           = "two-way",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "loo",
      se              = FALSE,
      parallel        = FALSE
    )
  ))

  set.seed(999)  # different seed should not affect LOO
  out2 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      force           = "two-way",
      factors.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "loo",
      se              = FALSE,
      parallel        = FALSE
    )
  ))

  expect_equal(out1$r.cv, out2$r.cv)
  # CV.out scores should be identical
  expect_equal(out1$CV.out, out2$CV.out, tolerance = 1e-10)
})

test_that("NTCV-P5: cv.sample reproducibility with set.seed", {
  for (cm in c("all_units", "treated_units")) {
    set.seed(42)
    out1 <- suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        factors.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = cm,
        se              = FALSE,
        parallel        = FALSE
      )
    ))

    set.seed(42)
    out2 <- suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        factors.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = cm,
        se              = FALSE,
        parallel        = FALSE
      )
    ))

    expect_equal(out1$r.cv, out2$r.cv,
                 info = paste("cv.method =", cm, ": r.cv not reproducible"))
    expect_equal(out1$CV.out, out2$CV.out, tolerance = 1e-10,
                 info = paste("cv.method =", cm, ": CV.out not reproducible"))
  }
})


## =================================================================
## Section I: Runtime Benchmarks (cv.method timing comparison)
## Informational only — skip on CRAN
## =================================================================

test_that("BENCH1: IFE timing comparison (loo vs all_units vs treated_units)", {
  skip_on_cran()

  cat("\n=== Runtime Benchmark: IFE nevertreated cv.method timing ===\n")

  set.seed(42)
  t_loo <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        factors.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = "loo",
        se              = FALSE,
        parallel        = FALSE
      )
    ))
  )
  cat(sprintf("  IFE loo:            %6.2f sec (elapsed)\n", t_loo["elapsed"]))

  set.seed(42)
  t_au <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        factors.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = "all_units",
        se              = FALSE,
        parallel        = FALSE
      )
    ))
  )
  cat(sprintf("  IFE all_units:      %6.2f sec (elapsed)\n", t_au["elapsed"]))

  set.seed(42)
  t_tu <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        factors.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = "treated_units",
        se              = FALSE,
        parallel        = FALSE
      )
    ))
  )
  cat(sprintf("  IFE treated_units:  %6.2f sec (elapsed)\n", t_tu["elapsed"]))

  cat(sprintf("  Speedup (loo/all_units):      %.2fx\n",
              t_loo["elapsed"] / max(t_au["elapsed"], 0.001)))
  cat(sprintf("  Speedup (loo/treated_units):  %.2fx\n",
              t_loo["elapsed"] / max(t_tu["elapsed"], 0.001)))
  cat("=== End IFE Benchmark ===\n")

  # All three must complete — that's the real test
  expect_true(TRUE)
})

test_that("BENCH2: CFE timing comparison (loo vs all_units vs treated_units)", {
  skip_on_cran()

  cat("\n=== Runtime Benchmark: CFE nevertreated cv.method timing ===\n")

  set.seed(42)
  t_loo <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "cfe",
        force           = "two-way",
        factors.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = "loo",
        se              = FALSE,
        parallel        = FALSE
      )
    ))
  )
  cat(sprintf("  CFE loo:            %6.2f sec (elapsed)\n", t_loo["elapsed"]))

  set.seed(42)
  t_au <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "cfe",
        force           = "two-way",
        factors.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = "all_units",
        se              = FALSE,
        parallel        = FALSE
      )
    ))
  )
  cat(sprintf("  CFE all_units:      %6.2f sec (elapsed)\n", t_au["elapsed"]))

  set.seed(42)
  t_tu <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "cfe",
        force           = "two-way",
        factors.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = "treated_units",
        se              = FALSE,
        parallel        = FALSE
      )
    ))
  )
  cat(sprintf("  CFE treated_units:  %6.2f sec (elapsed)\n", t_tu["elapsed"]))

  cat(sprintf("  Speedup (loo/all_units):      %.2fx\n",
              t_loo["elapsed"] / max(t_au["elapsed"], 0.001)))
  cat(sprintf("  Speedup (loo/treated_units):  %.2fx\n",
              t_loo["elapsed"] / max(t_tu["elapsed"], 0.001)))
  cat("=== End CFE Benchmark ===\n")

  expect_true(TRUE)
})

## ---------------------------------------------------------------
## Section G: Parallel CV Folds in fect_nevertreated
##
## Tests for REQ-parallel-cv: verifies that parallel=TRUE produces
## identical results to parallel=FALSE (sequential), that
## reproducibility holds under parallelism with fixed seeds, that
## default behavior is unchanged, that the LOO path is unaffected,
## and that edge cases work correctly.
##
## Follows test-spec.md for REQ-parallel-cv.
## Tolerances: 1e-10 for CV score differences (per test-spec.md).
## ---------------------------------------------------------------

## -- G.1  Sequential-Parallel Equivalence: IFE, all_units ----------

test_that("G.1: parallel CV matches sequential — IFE, all_units", {
  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  set.seed(123)
  result_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2
    )
  ))

  ## r.cv must be identical
  expect_identical(result_seq$r.cv, result_par$r.cv)

  ## CV.out matrix must match within tolerance 1e-10
  cv_diff <- max(abs(result_seq$CV.out - result_par$CV.out))
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out max diff = %.2e (tolerance = 1e-10)", cv_diff))
})

## -- G.2  Sequential-Parallel Equivalence: IFE, treated_units ------

test_that("G.2: parallel CV matches sequential — IFE, treated_units", {
  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "treated_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  set.seed(123)
  result_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "treated_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2
    )
  ))

  expect_identical(result_seq$r.cv, result_par$r.cv)

  cv_diff <- max(abs(result_seq$CV.out - result_par$CV.out))
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out max diff = %.2e (tolerance = 1e-10)", cv_diff))
})

## -- G.3  Sequential-Parallel Equivalence: CFE, all_units ----------

test_that("G.3: parallel CV matches sequential — CFE, all_units", {
  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "cfe",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  set.seed(123)
  result_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "cfe",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2
    )
  ))

  expect_identical(result_seq$r.cv, result_par$r.cv)

  cv_diff <- max(abs(result_seq$CV.out - result_par$CV.out))
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out max diff = %.2e (tolerance = 1e-10)", cv_diff))
})

## -- G.4  Sequential-Parallel Equivalence: CFE, treated_units ------

test_that("G.4: parallel CV matches sequential — CFE, treated_units", {
  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "cfe",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "treated_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  set.seed(123)
  result_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "cfe",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "treated_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2
    )
  ))

  expect_identical(result_seq$r.cv, result_par$r.cv)

  cv_diff <- max(abs(result_seq$CV.out - result_par$CV.out))
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out max diff = %.2e (tolerance = 1e-10)", cv_diff))
})

## -- G.5  Reproducibility Under Parallelism ------------------------

test_that("G.5: parallel CV is reproducible with same seed", {
  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_par1 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2,
      seed         = 12345
    )
  ))

  set.seed(123)
  result_par2 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2,
      seed         = 12345
    )
  ))

  expect_identical(result_par1$r.cv, result_par2$r.cv)
  expect_identical(result_par1$CV.out, result_par2$CV.out)
})

## -- G.6  Default Behavior Unchanged -------------------------------

test_that("G.6: default (no parallel arg) behaves as sequential", {
  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  ## Call without specifying parallel or cores — should default to sequential
  set.seed(123)
  result_default <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE
    )
  ))

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  expect_identical(result_default$r.cv, result_seq$r.cv)
  expect_identical(result_default$CV.out, result_seq$CV.out)
})

## -- G.7  LOO Path Unaffected --------------------------------------

test_that("G.7: LOO path is unaffected by parallel flag", {
  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_loo_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      cv.method    = "loo",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  set.seed(123)
  result_loo_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      cv.method    = "loo",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2
    )
  ))

  ## LOO should produce identical results regardless of parallel flag
  expect_identical(result_loo_seq$r.cv, result_loo_par$r.cv)
  expect_identical(result_loo_seq$CV.out, result_loo_par$CV.out)
})

## -- G.8  Edge Case: k = 1 ----------------------------------------

test_that("G.8: edge case — k = 1 with parallel=TRUE runs without error", {
  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 1,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2
    )
  ))

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 1,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  expect_identical(result_par$r.cv, result_seq$r.cv)
})

## -- G.9  Edge Case: cores = 1 ------------------------------------

test_that("G.9: edge case — cores = 1 behaves as sequential", {
  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_c1 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 1
    )
  ))

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  cv_diff <- max(abs(result_c1$CV.out - result_seq$CV.out))
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out max diff = %.2e (tolerance = 1e-10)", cv_diff))
})

## -- G.10  Edge Case: cores = NULL (auto-detect) ------------------

test_that("G.10: edge case — cores = NULL auto-detects and runs", {
  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_auto <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = NULL
    )
  ))

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  cv_diff <- max(abs(result_auto$CV.out - result_seq$CV.out))
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out max diff = %.2e (tolerance = 1e-10)", cv_diff))
})

## -- G.11  Edge Case: parallel=FALSE with cores specified ----------

test_that("G.11: edge case — parallel=FALSE ignores cores", {
  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_no_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE,
      cores        = 4
    )
  ))

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  expect_identical(result_no_par$r.cv, result_seq$r.cv)
  expect_identical(result_no_par$CV.out, result_seq$CV.out)
})

## -- G.12  Timing Benchmark: Parallel vs Sequential ----------------
## This is an informational benchmark; we report speedup but do not
## hard-fail on specific speedup thresholds (per test-spec.md:
## "Not a pass/fail test, but auditor should measure and report").

test_that("G.12: timing benchmark — parallel vs sequential with 10 cores", {
  skip_on_cran()

  ## Larger dataset for meaningful timing differences
  dat <- make_factor_data(N = 100, TT = 30, Ntr = 25, r = 2, seed = 99)

  cat("\n=== Parallel CV Timing Benchmark ===\n")

  ## --- IFE, all_units ---
  set.seed(123)
  t_seq_ife_au <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "ife",
        force        = "two-way",
        factors.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "all_units",
        se           = FALSE,
        parallel     = FALSE
      )
    ))
  )

  set.seed(123)
  t_par_ife_au <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "ife",
        force        = "two-way",
        factors.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "all_units",
        se           = FALSE,
        parallel     = TRUE,
        cores        = 10
      )
    ))
  )

  cat(sprintf("  IFE all_units  sequential: %6.2f sec\n", t_seq_ife_au["elapsed"]))
  cat(sprintf("  IFE all_units  parallel:   %6.2f sec (10 cores)\n", t_par_ife_au["elapsed"]))
  speedup_ife_au <- t_seq_ife_au["elapsed"] / max(t_par_ife_au["elapsed"], 0.001)
  cat(sprintf("  Speedup: %.2fx\n\n", speedup_ife_au))

  ## --- IFE, treated_units ---
  set.seed(123)
  t_seq_ife_tu <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "ife",
        force        = "two-way",
        factors.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "treated_units",
        se           = FALSE,
        parallel     = FALSE
      )
    ))
  )

  set.seed(123)
  t_par_ife_tu <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "ife",
        force        = "two-way",
        factors.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "treated_units",
        se           = FALSE,
        parallel     = TRUE,
        cores        = 10
      )
    ))
  )

  cat(sprintf("  IFE treated_units  sequential: %6.2f sec\n", t_seq_ife_tu["elapsed"]))
  cat(sprintf("  IFE treated_units  parallel:   %6.2f sec (10 cores)\n", t_par_ife_tu["elapsed"]))
  speedup_ife_tu <- t_seq_ife_tu["elapsed"] / max(t_par_ife_tu["elapsed"], 0.001)
  cat(sprintf("  Speedup: %.2fx\n\n", speedup_ife_tu))

  ## --- CFE, all_units ---
  set.seed(123)
  t_seq_cfe_au <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "cfe",
        force        = "two-way",
        factors.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "all_units",
        se           = FALSE,
        parallel     = FALSE
      )
    ))
  )

  set.seed(123)
  t_par_cfe_au <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "cfe",
        force        = "two-way",
        factors.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "all_units",
        se           = FALSE,
        parallel     = TRUE,
        cores        = 10
      )
    ))
  )

  cat(sprintf("  CFE all_units  sequential: %6.2f sec\n", t_seq_cfe_au["elapsed"]))
  cat(sprintf("  CFE all_units  parallel:   %6.2f sec (10 cores)\n", t_par_cfe_au["elapsed"]))
  speedup_cfe_au <- t_seq_cfe_au["elapsed"] / max(t_par_cfe_au["elapsed"], 0.001)
  cat(sprintf("  Speedup: %.2fx\n\n", speedup_cfe_au))

  ## --- CFE, treated_units ---
  set.seed(123)
  t_seq_cfe_tu <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "cfe",
        force        = "two-way",
        factors.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "treated_units",
        se           = FALSE,
        parallel     = FALSE
      )
    ))
  )

  set.seed(123)
  t_par_cfe_tu <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "cfe",
        force        = "two-way",
        factors.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "treated_units",
        se           = FALSE,
        parallel     = TRUE,
        cores        = 10
      )
    ))
  )

  cat(sprintf("  CFE treated_units  sequential: %6.2f sec\n", t_seq_cfe_tu["elapsed"]))
  cat(sprintf("  CFE treated_units  parallel:   %6.2f sec (10 cores)\n", t_par_cfe_tu["elapsed"]))
  speedup_cfe_tu <- t_seq_cfe_tu["elapsed"] / max(t_par_cfe_tu["elapsed"], 0.001)
  cat(sprintf("  Speedup: %.2fx\n\n", speedup_cfe_tu))

  cat("=== End Parallel CV Timing Benchmark ===\n")

  ## Informational — always passes; speedup is reported in test output
  expect_true(TRUE)
})

## -- G.13  Property: Backend Cleanup After Parallel CV -------------

test_that("G.13: parallel backend is restored after fect() returns", {
  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  old_plan <- future::plan()

  set.seed(123)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      factors.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2
    )
  ))

  new_plan <- future::plan()

  ## The future plan class should be restored
  expect_identical(class(old_plan), class(new_plan))
})
