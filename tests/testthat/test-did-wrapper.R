test_that("did_wrapper twfe runs and returns structure", {
  skip_if_not_installed("fixest")
  suppressWarnings(try(data("simdata", package = "fect"), silent = TRUE))
  expect_true(exists("simdata"))

  # create a binary ever-treated indicator for TWFE event study helper
  sim <- simdata
  sim$treat <- ave(sim$D, sim$id, FUN = function(z) as.numeric(mean(z, na.rm = TRUE) > 0))

  res <- did_wrapper(
    data = sim,
    Y = "Y",
    D = "D",
    X = c("X1", "X2"),
    index = c("id", "time"),
    method = "twfe",
    se = "default",
    parallel = FALSE
  )

  expect_s3_class(res, "did_wrapper")
  expect_true(is.data.frame(res$est.avg))
  expect_true(is.data.frame(res$est.att))
  expect_true(all(colnames(res$est.avg) %in% c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")))
})
