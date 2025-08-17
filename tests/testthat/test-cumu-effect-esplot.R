test_that("cumu/att.cumu/esplot run without error on fect output", {
  suppressWarnings(try(data("simdata", package = "fect"), silent = TRUE))
  expect_true(exists("simdata"))
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
  # cumu (att.cumu)
  # att.cumu requires uncertainty objects from keep.sims=TRUE
  out_keep <- fect::fect(
    Y ~ D + X1 + X2,
    data = simdata,
    index = c("id", "time"),
    method = "ife",
    r = 1,
    CV = FALSE,
    se = TRUE,
    nboots = 20,
    keep.sims = TRUE,
    parallel = FALSE
  )
  c <- att.cumu(out_keep, period = c(1, 3), plot = FALSE)
  expect_true(is.matrix(c) || is.data.frame(c))
  # esplot
  p <- esplot(out_keep$est.att, Estimate = "ATT")
  expect_s3_class(p, "ggplot")
})
