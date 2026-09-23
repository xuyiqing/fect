## Regression test: the LOO TOST must use the same pre-treatment window as the
## F test (periods passing the `proportion` cutoff). Before the fix, every row of
## pre.est.att entered the TOST, so a sparse early period with a wide S.E. could
## set the reported maximum p-value even though it was neither plotted nor part
## of the F test (found 2026-09-21 via the manual's MC LOO panel: 0.545 vs 0.0004).

test_that("LOO TOST p-value is computed over the proportion-filtered window", {
  set.seed(1)
  times <- -9:5                      # 10 pre-periods, 5 post-periods
  count <- c(3, 3, 60, 70, 80, 90, 100, 100, 100, 100, rep(100, 5))
  pre <- times[times <= 0]
  nboots <- 50
  ## pre-treatment averages: tight around zero, except two sparse early
  ## periods (3 treated units) with wide S.E.
  att <- rep(0.05, length(pre)); se <- rep(0.10, length(pre))
  att[1:2] <- c(-1.4, 1.3); se[1:2] <- c(1.5, 1.5)
  pre.est.att <- cbind(ATT = att, S.E. = se, CI.lower = att - 1.96 * se,
                       CI.upper = att + 1.96 * se, p.value = 0.5, count.on = count[1:10])
  rownames(pre.est.att) <- pre
  pre.att.boot <- matrix(rnorm(length(pre) * nboots, att, 0.10), length(pre), nboots)
  rownames(pre.att.boot) <- pre
  x <- list(time = times, count = count, sigma2.fect = 1,
            placeboTest = FALSE, carryoverTest = FALSE, loo = TRUE,
            pre.est.att = pre.est.att, pre.att.boot = pre.att.boot)
  out <- fect:::diagtest(x, proportion = 0.3, tost.threshold = 0.5)

  ## window = periods with count >= 0.3 * max(count) = 30, i.e. -7..0 (8 periods)
  expect_equal(out$df1, 8)
  ## every in-window period has |ATT| = 0.05, S.E. = 0.10 -> TOST p ~ 3e-6
  expect_lt(out$tost.equiv.p, 1e-4)
  ## the two sparse periods alone would give p ~ 0.5; they must not enter
  tost_sparse <- max(1 - pnorm((att[1:2] + 0.5) / se[1:2]),
                     1 - pnorm((0.5 - att[1:2]) / se[1:2]))
  expect_gt(tost_sparse, 0.4)
})
