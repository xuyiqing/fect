## ---------------------------------------------------------------
## Regression test: Y.ct.full populated at control positions on GSC.
##
## Pre-fix (HANDOFF 2026-04-25): On the GSC path
## (method = "ife" + time.component.from = "nevertreated"), Y.ct.full
## at id_co positions was set via `Y.co - residuals`, which leaves NA
## at any position where Y.co is NA. This blocks user-space rolling CV
## from scoring MSPE at masked control observations -- the channel
## that R/cv-rolling.R uses to evaluate held-out positions.
##
## Fix (Stage 2, 2026-04-25): After Y.ct.full <- Y.ct on the shared
## post-bifurcation path of fect_nevertreated.R, on the IFE method
## branch (the GSC estimator) overwrite control columns with the
## model-implied factor product:
##     Y.ct.full[, id_co] <- F * t(lambda_co)
## sourced directly from est.co.best$factor and est.co.best$lambda.
##
## Tests:
##   1. fit$Y.ct.full[, fit$co] matches F * t(lambda) element-wise.
##   2. ATT.avg matches the unmodified-dev baseline (4.639593 on
##      simgsynth, set.seed(11), r=2, CV=FALSE, force="two-way").
## ---------------------------------------------------------------

test_that("GSC populates Y.ct.full[, co] with F * t(lambda)", {

  skip_on_cran()

  if (!exists("simgsynth", inherits = TRUE)) {
    suppressWarnings(try(utils::data("simgsynth", package = "fect"),
                         silent = TRUE))
  }
  skip_if_not(exists("simgsynth", inherits = TRUE),
              "simgsynth dataset not available")

  set.seed(11)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simgsynth,
      index     = c("id", "time"),
      method    = "ife",
      time.component.from = "nevertreated",
      r         = 2,
      CV        = FALSE,
      force     = "two-way",
      se        = FALSE
    )
  ))

  expect_true(!is.null(fit$factor))
  expect_true(!is.null(fit$lambda.co))
  expect_true(!is.null(fit$Y.ct.full))
  expect_true(!is.null(fit$co))

  F.hat <- as.matrix(fit$factor)
  L.co  <- as.matrix(fit$lambda.co)
  expected <- F.hat %*% t(L.co)

  observed <- fit$Y.ct.full[, fit$co, drop = FALSE]

  expect_equal(dim(observed), dim(expected))
  expect_equal(observed, expected, tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("GSC ATT.avg unchanged after Stage 2 patch (simgsynth anchor)", {

  skip_on_cran()

  if (!exists("simgsynth", inherits = TRUE)) {
    suppressWarnings(try(utils::data("simgsynth", package = "fect"),
                         silent = TRUE))
  }
  skip_if_not(exists("simgsynth", inherits = TRUE),
              "simgsynth dataset not available")

  set.seed(11)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simgsynth,
      index     = c("id", "time"),
      method    = "ife",
      time.component.from = "nevertreated",
      r         = 2,
      CV        = FALSE,
      force     = "two-way",
      se        = FALSE
    )
  ))

  ## Anchor: 4.639593 on unmodified dev, set.seed(11), simgsynth.
  expect_equal(fit$att.avg, 4.639593, tolerance = 1e-5)
})
