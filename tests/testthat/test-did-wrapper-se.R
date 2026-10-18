## did_wrapper standard errors (2.4.8): the default call, the cluster
## bootstrap draw, and the jackknife scaling.

test_that("did_wrapper runs with method and se left at their defaults", {
  skip_on_cran()
  skip_if_not_installed("fixest")
  data("simdata", package = "fect", envir = environment())
  res <- suppressMessages(did_wrapper(simdata, "Y", "D", index = c("id", "time")))
  expect_s3_class(res, "did_wrapper")
  expect_true(is.finite(res$est.avg$ATT.avg))
})

test_that("a cluster-bootstrap draw keeps repeated units as separate units", {
  d <- data.frame(id = rep(c("a", "b", "c", "d", "e"), each = 3), t = rep(1:3, 5), y = 1:15)
  ids <- unique(d$id)
  rows <- split(seq_len(nrow(d)), factor(d$id, levels = ids))
  set.seed(7)
  for (b in 1:20) {
    d_b <- fect:::.did_cluster_boot_draw(d, "id", rows)
    expect_equal(nrow(d_b), nrow(d))                  # every draw has 5 units x 3 periods
    expect_equal(length(unique(d_b$id)), length(ids)) # copies get their own IDs
    expect_true(all(table(d_b$id) == 3))
  }
  ## with 5 units, some draw repeats a unit: its rows appear more than once
  set.seed(7)
  dup <- vapply(1:20, function(b) {
    d_b <- fect:::.did_cluster_boot_draw(d, "id", rows)
    anyDuplicated(d_b$y) > 0
  }, logical(1))
  expect_true(any(dup))
})

test_that("did_wrapper jackknife SE uses the (n - 1) / n scaling", {
  skip_on_cran()
  skip_if_not_installed("fixest")
  data("simdata", package = "fect", envir = environment())
  sim <- simdata[simdata$id %in% sort(unique(simdata$id))[1:60], ]
  res <- suppressMessages(suppressWarnings(
    did_wrapper(sim, "Y", "D", index = c("id", "time"), method = "twfe",
                se = "jackknife", parallel = FALSE)
  ))
  ids <- unique(sim$id)
  loo <- vapply(ids, function(i) {
    fit <- fixest::feols(Y ~ D | id + time, data = sim[sim$id != i, ], notes = FALSE)
    unname(stats::coef(fit)["D"])
  }, numeric(1))
  n <- length(loo)
  se_jack <- sqrt((n - 1) / n * sum((loo - mean(loo))^2))
  expect_equal(res$est.avg$S.E., se_jack, tolerance = 1e-8)
  expect_equal(res$est.avg$CI.lower, res$est.avg$ATT.avg - 1.96 * se_jack, tolerance = 1e-8)
  ## the jackknife SE is close to the cluster-robust SE, not sqrt(n) smaller
  res0 <- suppressMessages(did_wrapper(sim, "Y", "D", index = c("id", "time"),
                                       method = "twfe", se = "default", parallel = FALSE))
  expect_gt(res$est.avg$S.E. / res0$est.avg$S.E., 0.7)
  expect_lt(res$est.avg$S.E. / res0$est.avg$S.E., 1.4)
})
