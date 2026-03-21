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

test_that("plot.fect executes for heterogeneous (hte) plot without error", {
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
  p <- plot.fect(out, type = "hte", covariate = "X1", show.count = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot.fect executes for calendar plot with relative.time=TRUE without error", {
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
  p <- plot.fect(out, type = "calendar", relative.time = TRUE, show.count = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plot.fect executes for heterogeneous pretreatment plot with num.pretreatment without error", {
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
  p <- plot.fect(out, type = "hte", covariate = "X1", pretreatment = TRUE, num.pretreatment = 3, show.count = FALSE)
  expect_s3_class(p, "ggplot")
})