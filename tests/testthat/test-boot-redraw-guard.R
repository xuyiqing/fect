## v2.5.1: the bootstrap redraw loops stop after a fixed number of tries
## instead of running forever on a sparse panel.

test_that("the redraw counter stops with a classed error after max.tries", {
  tick <- fect:::.fect_redraw_guard(max.tries = 3L)
  expect_silent(for (i in 1:3) tick())
  expect_error(tick(), class = "fect_redraw_error")
})

test_that("bootstrap on a panel no resample can cover ends with an error", {
  skip_on_cran()
  ## 36 control units each observed in 5 periods, treated units treated
  ## from period 21: each of periods 21-200 has one untreated observation,
  ## so a resample covers every period with probability near 1e-7.
  rows <- list()
  for (u in 1:40) rows[[u]] <- data.frame(id = u, time = ((u - 1) * 5 + 1):(u * 5), D = 0L)
  for (u in 41:50) rows[[u]] <- data.frame(id = u, time = 1:200, D = as.integer(1:200 > 20))
  d <- do.call(rbind, rows)
  set.seed(1)
  d$Y <- rnorm(nrow(d)) + d$D
  t0 <- Sys.time()
  expect_error(
    suppressWarnings(suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"),
      method = "fe", se = TRUE, nboots = 20, parallel = FALSE, min.T0 = 1))),
    class = "fect_redraw_error"
  )
  expect_lt(as.numeric(Sys.time() - t0, units = "secs"), 60)
})
