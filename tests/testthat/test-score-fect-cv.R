## ---------------------------------------------------------------
## Tests for fect_cv: regression coverage (S2) + cv.method extension
## in fect_cv (Section B, NEW for Phase 2).
##
## Originally part of test-score-unify.R; split out 2026-05-03 for
## progress visibility. Shared fixtures live in helper-score-unify.R.
## ---------------------------------------------------------------

## =================================================================
## S2: fect_cv regression test (before vs after refactor)
## =================================================================

test_that("S2.1: IFE method CV - r.cv and CV.out snapshot", {

  skip_on_cran()
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

  skip_on_cran()
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

  skip_on_cran()
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
## Section B: cv.method in fect_cv (NEW for Phase 2)
## =================================================================

test_that("CV1: cv.method='all_units' selects r.cv", {

  skip_on_cran()
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

  skip_on_cran()
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

  skip_on_cran()
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


