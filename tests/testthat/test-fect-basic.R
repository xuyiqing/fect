test_that("fect formula basic runs on simdata and returns expected slots", {
  suppressWarnings(try(data("simdata", package = "fect"), silent = TRUE))
  expect_true(exists("simdata"))

  set.seed(1)
  out <- fect::fect(
    Y ~ D + X1 + X2,
    data = simdata,
    index = c("id", "time"),
    method = "ife",
    r = 2,
    CV = FALSE,
    se = FALSE,
    parallel = FALSE
  )

  expect_s3_class(out, "fect")
  expect_true(is.matrix(out$eff))
  expect_true(is.numeric(out$att.avg))
  expect_true(length(out$time) >= 1)
})

test_that("fect returns est.att and att.vcov when se=TRUE", {
  suppressWarnings(try(data("simdata", package = "fect"), silent = TRUE))
  expect_true(exists("simdata"))

  set.seed(2)
  out <- fect::fect(
    Y ~ D + X1 + X2,
    data = simdata,
    index = c("id", "time"),
    method = "ife",
    r = 1,
    CV = FALSE,
    se = TRUE,
    nboots = 30,
    parallel = FALSE
  )

  expect_true(!is.null(out$est.att))
  expect_true(!is.null(out$att.vcov))
  expect_true(all(colnames(out$est.att) %in% c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value", "count")))
})
