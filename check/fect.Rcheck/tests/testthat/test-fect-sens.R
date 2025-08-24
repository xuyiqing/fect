test_that("fect_sens attaches sensitivity results when inputs present", {
  skip_if_not_installed("HonestDiDFEct")
  suppressWarnings(try(data("simdata", package = "fect"), silent = TRUE))
  expect_true(exists("simdata"))

  set.seed(3)
  out <- fect::fect(
    Y ~ D + X1 + X2,
    data = simdata,
    index = c("id", "time"),
    method = "ife",
    r = 1,
    CV = FALSE,
    se = TRUE,
    nboots = 20,
    parallel = FALSE
  )

  # ensure names and shapes needed by fect_sens exist
  expect_true(!is.null(out$est.att))
  expect_true(!is.null(out$att.vcov))

  # Ensure placebo.period exists to avoid length-0 error inside fect_sens
  if (is.null(out$placebo.period)) {
    out$placebo.period <- c(-3, -1)
  }
  # restrict periods so numPre+numPost matches available row count
  out$placebo.period <- c(-3, -1)
  out2 <- fect_sens(out, post.periods = 1:10, parallel = FALSE)
  expect_true(!is.null(out2$sensitivity.rm) || !is.null(out2$sensitivity.smooth))
})
