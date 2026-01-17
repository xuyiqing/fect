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

  res <- plot.fect(out, type = "gap", show.count = FALSE, return.data = TRUE)
  expect_true(is.list(res))
  expect_s3_class(res, "fect_plot_return")
  expect_s3_class(res$p, "ggplot")
  expect_true(is.list(res$data))
  expect_true(is.data.frame(res$data$estimate))
  expect_true(all(c("Period", "ATT", "plot_type") %in% names(res$data$estimate)))
  expect_true(all(res$data$estimate$plot_type == "gap"))
})
