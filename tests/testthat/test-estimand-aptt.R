## ---------------------------------------------------------------
## Tests for estimand(fit, "aptt", "event.time").
##
## Reference: Chen & Roth (2024 QJE).
## APTT_t = mean_{(t,i): D=1} (Y - Y0_hat) / mean_{(t,i): D=1} Y0_hat
##
## Validates against a manual computation from raw fit slots — this
## is the same recipe currently in vignettes/02-fect.Rmd.
## ---------------------------------------------------------------

suppressWarnings(data("simdata", package = "fect"))


.fit_canonical <- function(nboots = 50, keep.sims = TRUE) {
  set.seed(42)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = simdata, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, nboots = nboots, parallel = FALSE,
      keep.sims = keep.sims
    )
  ))
}


## -- AP.1  Schema --------------------------------------------------

test_that("AP.1: estimand(fit, 'aptt', 'event.time') has the documented schema", {

  skip_on_cran()

  fit <- .fit_canonical()
  est <- fect::estimand(fit, "aptt", "event.time")

  expect_s3_class(est, "data.frame")
  expected_cols <- c("event.time", "estimate", "se",
                     "ci.lo", "ci.hi", "n_cells", "vartype")
  expect_setequal(names(est), expected_cols)

  expect_true(is.numeric(est$estimate))
  expect_true(all(!is.na(est$event.time)))
})


## -- AP.2  Manual recipe match -----------------------------------------

test_that("AP.2: estimate matches the manual APTT recipe from the vignette", {

  skip_on_cran()

  fit <- .fit_canonical()
  est <- fect::estimand(fit, "aptt", "event.time")

  ## Manual APTT computation: per event time, mean(eff)/mean(Y - eff).
  Don <- !is.na(fit$D.dat) & fit$D.dat == 1 & !is.na(fit$T.on)
  ets <- sort(unique(fit$T.on[Don]))

  manual <- vapply(ets, function(et) {
    m <- Don & fit$T.on == et
    eff_t <- fit$eff[m]
    Y_t   <- fit$Y.dat[m]
    Y0_t  <- Y_t - eff_t
    mean(eff_t, na.rm = TRUE) / mean(Y0_t, na.rm = TRUE)
  }, numeric(1))

  est_sorted <- est[order(est$event.time), ]
  expect_equal(est_sorted$estimate, manual, tolerance = 1e-10)
})


## -- AP.3  keep.sims = FALSE: error path ------------------------------

test_that("AP.3: APTT without keep.sims errors with the locked wording", {

  skip_on_cran()

  fit_no_keep <- .fit_canonical(keep.sims = FALSE)
  expect_error(
    fect::estimand(fit_no_keep, "aptt", "event.time"),
    "keep\\.sims = TRUE",
    fixed = FALSE
  )
})


## -- AP.4  vartype = "none" works without keep.sims --------------------

test_that("AP.4: vartype = 'none' returns point estimate only without keep.sims", {

  skip_on_cran()

  fit_no_keep <- .fit_canonical(keep.sims = FALSE)
  est <- fect::estimand(fit_no_keep, "aptt", "event.time", vartype = "none")
  expect_true(nrow(est) > 0)
  expect_true(all(is.na(est$se)))
  expect_true(all(is.na(est$ci.lo)))
  expect_true(all(est$vartype == "none"))
})
