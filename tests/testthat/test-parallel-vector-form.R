## ---------------------------------------------------------------
## Regression test: parallel = c("cv", "boot") must not error.
##
## Phase 1-3 updated the CV dispatch (cv.R, fect_nevertreated.R) to
## handle the vector form via the do_parallel_cv / do_parallel_boot
## flags in default.R. Several downstream paths (fect_boot, fit_test,
## permutation, fect_sens, did_wrapper) kept legacy `parallel == TRUE`
## or `if (parallel)` checks that error with
##   "the condition has length > 1"
## when the user passes a length > 1 character vector.
##
## Caught by an Alsaadi (2025) real-workload smoke test on 2026-04-22.
## ---------------------------------------------------------------

suppressWarnings(data("simdata", package = "fect"))

test_that("parallel = c('cv','boot') does not crash fect_boot", {

  skip_on_cran()

  set.seed(42)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:1,
      CV        = TRUE,
      k         = 2,
      cv.method = "all_units",
      time.component.from = "nevertreated",
      se        = TRUE,
      vartype   = "bootstrap",
      nboots    = 5,
      parallel  = c("cv", "boot"),
      cores     = 2
    )
  ))

  expect_true(!is.null(fit$att.avg))
  expect_true(!is.null(fit$est.att))
})
