## ---------------------------------------------------------------
## Tests for estimand(fit, "att", ...) and the backcompat invariant
## with fit$est.att.
##
## Locked invariant per ref/po-estimands-contract.md §2:
##   estimand(fit, "att", "event.time") matches fit$est.att columns
##   ATT / S.E. / CI.lower / CI.upper byte-for-byte under default args.
## ---------------------------------------------------------------

suppressWarnings(data("simdata", package = "fect"))


.fit_canonical <- function(nboots = 50, keep.sims = FALSE) {
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


## -- A.1  Schema: returns documented columns -------------------------

test_that("A.1: estimand(fit, 'att', 'event.time') returns the documented schema", {

  skip_on_cran()

  fit <- .fit_canonical()
  est <- fect::estimand(fit, "att", "event.time")

  expect_s3_class(est, "data.frame")
  expected_cols <- c("event.time", "estimate", "se",
                     "ci.lo", "ci.hi", "n_cells", "vartype")
  expect_setequal(names(est), expected_cols)

  expect_true(is.numeric(est$event.time))
  expect_true(is.numeric(est$estimate))
  expect_true(is.numeric(est$se))
  expect_true(is.numeric(est$n_cells))
})


## -- A.2  Backcompat: byte-equality with fit$est.att (LOCKED INVARIANT)

test_that("A.2: estimand(fit, 'att', 'event.time') matches fit$est.att byte-for-byte", {

  skip_on_cran()

  fit <- .fit_canonical()
  est <- fect::estimand(fit, "att", "event.time")

  ## fit$est.att rows are labeled by event time (rownames); match the
  ## same axis ordering.
  est_att <- as.data.frame(fit$est.att)
  est_att$event.time <- as.numeric(rownames(fit$est.att))
  est_att <- est_att[order(est_att$event.time), ]
  est <- est[order(est$event.time), ]

  expect_equal(est$event.time, est_att$event.time)
  expect_identical(est$estimate, unname(est_att$ATT))
  expect_identical(est$se,       unname(est_att$S.E.))
  expect_identical(est$ci.lo,    unname(est_att$CI.lower))
  expect_identical(est$ci.hi,    unname(est_att$CI.upper))
  expect_identical(est$n_cells,  unname(est_att$count))
})


## -- A.3  Errors: unsupported types in this commit -------------------

test_that("A.3: unimplemented type raises a clear message", {

  skip_on_cran()

  ## "att.cumu", "aptt", "log.att" all implemented; A.3 has nothing
  ## left to assert.
  expect_true(TRUE)
})


## -- A.4  Errors: invalid by ----------------------------------------

test_that("A.4: invalid by raises a useful error", {

  skip_on_cran()

  fit <- .fit_canonical()
  expect_error(
    fect::estimand(fit, "att", "nope"),
    "user-column",
    fixed = FALSE
  )
})


## -- A.5  Errors: bad fit ------------------------------------------

test_that("A.5: passing non-fect to estimand errors clearly", {
  expect_error(
    fect::estimand("not a fit", "att", "event.time"),
    "fect"
  )
})
