test_that("plot.fect executes for gap plot without error", {
  suppressWarnings(try(data("simdata", package = "fect"), silent = TRUE))
  expect_true(exists("simdata"))
  set.seed(4)
  out <- fect::fect(
    Y ~ D + X1 + X2,
    data = simdata,
    index = c("id", "time"),
    method = "ife",
    r = 1,
    CV = FALSE,
    se = FALSE,
    parallel = FALSE
  )
  # Just ensure the plotting call produces a ggplot object
  p <- plot.fect(out, type = "gap", show.count = FALSE)
  expect_s3_class(p, "ggplot")
})
