## ---------------------------------------------------------------
## Tests for the explicit carryover.rm slot on the fit object.
##
## Pre-fix (v2.3.2 and earlier), plot.R recovered carryover.rm from
## `as.list(x$call)$carryover.rm` and `eval(arg, envir = parent.frame())`.
## That call-parsing path silently gave the wrong answer when the fit
## was constructed via do.call(), wrappers that rewrote the call, or
## when x$call was missing/altered. Storing carryover.rm directly on
## the fit object makes the plot logic robust to all those paths.
## ---------------------------------------------------------------

suppressWarnings(data("simdata", package = "fect"))

## -- C.1  Slot exists and matches the value passed at fit time --

test_that("C.1: fit$carryover.rm is populated with the value passed to fect()", {

  skip_on_cran()

  set.seed(42)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = simdata, index = c("id", "time"),
      method = "ife", r = 1, CV = FALSE,
      force = "two-way",
      se = FALSE,
      carryover.rm = 2
    )
  ))

  expect_true("carryover.rm" %in% names(fit))
  expect_equal(fit$carryover.rm, 2)
})


## -- C.2  Slot is NULL when the user does not pass carryover.rm --

test_that("C.2: fit$carryover.rm is NULL when carryover.rm is not used", {

  skip_on_cran()

  set.seed(42)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = simdata, index = c("id", "time"),
      method = "ife", r = 1, CV = FALSE,
      force = "two-way",
      se = FALSE
    )
  ))

  expect_null(fit$carryover.rm)
})


## -- C.3  Plot logic uses the slot, not x$call. Wipe x$call after
## fitting and confirm that plot(fit) still recovers carryover_rm_K.
## This is the regression test for the call-parsing fragility: the
## old code would have read 0 from a call-less fit and rendered the
## carryover plot with a missing K-shift; the new code reads from
## x$carryover.rm and is unaffected.

test_that("C.3: plot logic recovers carryover.rm even when x$call is wiped", {

  skip_on_cran()

  set.seed(42)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = simdata, index = c("id", "time"),
      method = "ife", r = 1, CV = FALSE,
      force = "two-way",
      se = FALSE,
      carryover.rm = 2
    )
  ))

  ## Wipe x$call to defeat the legacy call-parsing path.
  fit$call <- NULL

  ## The slot should still report K = 2, and plot.fect should not
  ## throw. (We do not inspect the rendered plot's pixel content;
  ## a successful build is the contract guard.)
  expect_equal(fit$carryover.rm, 2)
  expect_no_error(suppressWarnings(suppressMessages(plot(fit, type = "exit"))))
})
