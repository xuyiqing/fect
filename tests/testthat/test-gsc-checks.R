## fect 2.5.0: the gsynth-revisit steps (R/gsc-checks.R, R/gsc-rank.R,
## R/ife-boot.R). Fast unit tests; the comparisons with the paper's
## replication code are in test-gsc-reference.R.

data(simgsynth, package = "fect")

## a treated unit far outside the controls' loadings: unit 101's outcome
## gets 20 times the first factor
.far_data <- function() {
  d <- simgsynth
  k <- d$id == 101
  d$Y[k] <- d$Y[k] + 20 * d$F1[k]
  d
}

test_that("least squares on the simplex: the active set solves the dense problem", {
  set.seed(3)
  y <- rnorm(12)
  X <- matrix(rnorm(12 * 260), 12, 260)
  w1 <- fect:::.gsc_simplex_ls(y, X, dense.max = 0L)       # active set
  w2 <- fect:::.gsc_simplex_ls(y, X, dense.max = 1000L)    # dense quadprog
  expect_equal(sum(w1), 1, tolerance = 1e-10)
  expect_true(all(w1 >= -1e-12))
  f <- function(w) sum((y - X %*% w)^2)
  expect_equal(f(w1), f(w2), tolerance = 1e-7)
  ## optimality: no column improves on the multiplier
  g <- drop(crossprod(X, X %*% w1 - y))
  expect_gte(min(g - sum(w1 * g)), -1e-8)
})

test_that("overlap_check(): tables, weights, flags and the caller's RNG", {
  set.seed(42); a <- runif(1)
  set.seed(42)
  oc <- overlap_check(Y ~ D, data = simgsynth, index = c("id", "time"))
  b <- runif(1)
  expect_identical(a, b)                                   # RNG state restored
  expect_s3_class(oc, "fect_overlap")
  u <- oc$units
  expect_identical(nrow(u), 5L)
  expect_equal(sum(u$weight), 1)
  expect_true(all(u$status == "ok"))
  expect_true(all(u$pct >= 0 & u$pct <= 1))
  expect_identical(u$flag, u$error > u$yardstick)
  expect_identical(oc$summary$n.flagged, sum(u$flag))
  expect_identical(u$n.donors, rep(45L, 5))
  expect_identical(oc$level, 0.90)
  expect_true(all(u$onset == 21))
  expect_output(print(oc), "Overlap check without a rank")
  expect_s3_class(plot(oc), "ggplot")
})

test_that("overlap_check(): a level shift of a treated unit changes nothing", {
  d <- simgsynth
  d$Y[d$id == 101] <- d$Y[d$id == 101] + 100
  a <- overlap_check(Y ~ D, data = simgsynth, index = c("id", "time"))$units
  b <- overlap_check(Y ~ D, data = d, index = c("id", "time"))$units
  expect_equal(a$error, b$error, tolerance = 1e-10)
  expect_identical(a$flag, b$flag)
})

test_that("overlap_check(): a unit off the donors' span is flagged", {
  oc <- overlap_check(Y ~ D, data = .far_data(), index = c("id", "time"))
  expect_true(oc$units$flag[oc$units$unit == "101"])
  expect_identical(oc$units$pct[oc$units$unit == "101"], 1)
})

test_that("overlap_check(): not assessable, and input checks", {
  d <- simgsynth[simgsynth$time > 18, ]                    # two pretreatment periods
  oc <- overlap_check(Y ~ D, data = d, index = c("id", "time"))
  expect_true(all(oc$units$status == "not assessable"))
  expect_true(all(oc$units$reason == "too few periods"))
  expect_true(all(is.na(oc$units$flag)))
  expect_error(overlap_check(Y ~ D, data = simgsynth, index = c("id", "time"), level = 1), "level")
  expect_error(overlap_check(Y ~ D, data = simgsynth, index = c("id", "nope")), "index")
  expect_error(overlap_check(Y ~ Z, data = simgsynth, index = c("id", "time")), "not in the data")
  d2 <- simgsynth
  d2$D[d2$id == 101 & d2$time == 30] <- 0                  # switches off
  expect_error(overlap_check(Y ~ D, data = d2, index = c("id", "time")), "absorbing")
})

test_that("fect_rank(): picks, table and the forward window", {
  rk <- fect_rank(Y ~ D, data = simgsynth, index = c("id", "time"), r = 0:3, keep.fits = TRUE)
  expect_s3_class(rk, "fect_rank")
  expect_identical(rk$table$r, 0:3)
  expect_true(all(c("fwd", "loo", "ATT", "extrapolation") %in% names(rk$table)))
  expect_true(rk$r.fwd %in% 0:3)
  expect_true(rk$r.loo %in% 0:3)
  expect_identical(rk$r.fwd, rk$table$r[which.min(rk$table$fwd)])
  expect_identical(rk$r.loo, rk$table$r[which.min(rk$table$loo)])
  expect_true(rk$window$meets)                             # 20 pretreatment periods
  expect_identical(rk$n.fwd, 5L)
  ## the ATT in the table is the fit's at that rank
  f2 <- fect(Y ~ D, data = simgsynth, index = c("id", "time"), method = "gsynth", r = 2, se = FALSE)
  expect_equal(rk$table$ATT[rk$table$r == 2], unname(f2$att.avg))
  expect_output(print(rk), "Forward rule")
  expect_error(fect_rank(Y ~ D, data = simgsynth, index = c("id", "time"), method = "ife"), "sets")
  expect_error(fect_rank(Y ~ D, data = simgsynth, index = c("id", "time"), r = c(1, 2, 4)), "without gaps")
})

test_that("fect_rank(): a short window is flagged, and units with too few periods are left out", {
  d <- simgsynth[simgsynth$time > 8, ]                     # 12 pretreatment periods
  rk <- fect_rank(Y ~ D, data = d, index = c("id", "time"), r = 0:5, min.T0 = 5)
  expect_false(rk$window$meets)                            # 8 forward predictions need 15
  expect_identical(rk$rbar.fwd, 5L)
  d <- simgsynth[simgsynth$time > 14, ]                    # 6 pretreatment periods
  rk <- fect_rank(Y ~ D, data = d, index = c("id", "time"), r = 0:5, min.T0 = 5)
  expect_identical(rk$rbar.fwd, 3L)                        # 6 - (3 + 1) - 1 = 1 forward task
  expect_true(all(is.na(rk$table$fwd[rk$table$r > 3])))
})

test_that("fect(cv.method = \"forward\") fits at the forward rule's rank", {
  rk <- fect_rank(Y ~ D, data = simgsynth, index = c("id", "time"), r = 0:4)
  f <- suppressMessages(fect(Y ~ D, data = simgsynth, index = c("id", "time"), method = "gsynth",
                             CV = TRUE, r = c(0, 4), cv.method = "forward"))
  expect_identical(as.integer(f$r.cv), rk$r.fwd)
  expect_identical(f$rank.forward$table$fwd, rk$table$fwd)
  g <- fect(Y ~ D, data = simgsynth, index = c("id", "time"), method = "gsynth", r = rk$r.fwd, se = FALSE)
  expect_equal(f$att.avg, g$att.avg)
  expect_identical(f$cv.method, "forward")
  expect_error(fect(Y ~ D, data = simgsynth, index = c("id", "time"), method = "ife",
                    CV = TRUE, r = c(0, 3), cv.method = "forward"), "GSC")
  expect_error(fect(Y ~ D, data = simgsynth, index = c("id", "time"), method = "gsynth",
                    CV = FALSE, r = 2, cv.method = "forward"), "CV = TRUE")
})

test_that("loading_check(): hull, distances, flags and extrapolation", {
  fit <- fect(Y ~ D, data = simgsynth, index = c("id", "time"), method = "gsynth", r = 2, se = FALSE)
  lc <- loading_check(fit)
  expect_s3_class(lc, "fect_loading")
  t <- lc$table
  expect_identical(nrow(t), 5L)
  expect_true(all(t$d2[t$in.hull] == 0))
  expect_true(all(t$d2[!t$in.hull] > 0))
  expect_true(all(!t$flagged[t$in.hull]))
  expect_identical(t$flagged, !t$in.hull & t$d2 > qchisq(0.95, 2))
  expect_identical(nrow(lc$Lco), 45L)
  ## the extrapolation share as defined in the paper's summary code
  post <- fit$D.dat[, fit$tr] == 1
  ct <- fit$Y.ct[, fit$tr][post]
  rg <- range(fit$Y.dat[fit$I.dat == 1])
  expect_identical(lc$extrapolation$share, mean(ct < rg[1] | ct > rg[2]))
  expect_output(print(lc), "Loading check at r = 2")
  expect_s3_class(plot(lc), "ggplot")
})

test_that("loading_check(): inputs", {
  f0 <- fect(Y ~ D, data = simgsynth, index = c("id", "time"), method = "gsynth", r = 0, se = FALSE)
  expect_error(loading_check(f0), "no factors")
  fi <- fect(Y ~ D, data = simgsynth, index = c("id", "time"), method = "ife", r = 2, se = FALSE)
  expect_error(loading_check(fi), "GSC")
  expect_error(loading_check(list()), "fit of fect")
})

test_that("loading_check(trim = TRUE) drops the flagged unit, chooses the rank again and refits", {
  skip_on_cran()
  d <- .far_data()
  fit <- fect(Y ~ D, data = d, index = c("id", "time"), method = "gsynth", r = 2, se = FALSE)
  lc <- loading_check(fit, trim = TRUE)
  expect_true("101" %in% lc$trim$out)
  expect_identical(lc$trim$status, "ok")
  expect_identical(lc$trim$est$n.tr, 5L - length(lc$trim$out))
  ## the trimmed fit is the fit at the reselected rank on the kept units
  kept <- d[!as.character(d$id) %in% lc$trim$out, ]
  rk <- fect_rank(Y ~ D, data = kept, index = c("id", "time"), r = 0:5)
  expect_identical(lc$trim$r, rk$r.fwd)
  g <- fect(Y ~ D, data = kept, index = c("id", "time"), method = "gsynth", r = rk$r.fwd, se = FALSE)
  expect_equal(lc$trim$est$ATT, unname(g$att.avg))
  ## data given explicitly
  lc2 <- loading_check(fit, trim = TRUE, data = d)
  expect_identical(lc2$trim$est, lc$trim$est)
})

test_that("IFE-EM bootstrap: runs, records its settings, and is reproducible", {
  skip_on_cran()
  run <- function(...) suppressWarnings(suppressMessages(
    fect(Y ~ D + X1, data = simgsynth, index = c("id", "time"), method = "ife", r = 2, CV = FALSE,
         se = TRUE, vartype = "parametric", nboots = 20, seed = 7, keep.sims = TRUE, ...)))
  f1 <- run(parallel = FALSE)
  ib <- f1$ife.boot
  expect_identical(ib$r, 2L)
  expect_identical(ib$r.noise, 1L)                         # 45 controls: one factor fewer
  expect_identical(ib$G, 30L)
  expect_identical(length(ib$groups), 45L)
  expect_identical(sort(unique(ib$groups)), 1:30)
  expect_identical(dim(ib$noise), c(30L, 45L))
  expect_identical(dim(ib$draws), c(20L, 50L))
  expect_true(all(is.finite(f1$att.avg.boot)))
  expect_true(is.finite(f1$est.avg[1, "S.E."]))
  expect_identical(f1$vartype, "parametric")
  ## the draws center on zero (no effect in the bootstrap panel), not on the estimate
  expect_lt(abs(mean(f1$att.avg.boot)), abs(f1$att.avg))
  f2 <- run(parallel = FALSE)
  expect_identical(f2$att.avg.boot, f1$att.avg.boot)
  f3 <- run(parallel = TRUE, cores = 2)
  expect_identical(f3$ife.boot$draws, f1$ife.boot$draws)
  expect_equal(f3$att.avg.boot, f1$att.avg.boot, tolerance = 1e-10)
})

test_that("IFE-EM bootstrap: noise vectors are leave-group-out residuals", {
  skip_on_cran()
  set.seed(1)
  T <- 12; N <- 14; r <- 1
  F <- matrix(rnorm(T), T, 1); L <- matrix(rnorm(N), N, 1)
  Y <- F %*% t(L) + outer(rnorm(T), rep(1, N)) + outer(rep(1, T), rnorm(N)) + matrix(rnorm(T * N), T, N)
  D <- matrix(0, T, N); D[9:T, 1:2] <- 1
  I <- matrix(1, T, N); II <- (D == 0) * 1
  T.on <- matrix(NA, T, N)
  co <- 3:N
  set.seed(5)
  out <- fect:::.ife_noise_group(Y, NULL, D, I, II, T.on, co, r.noise = 0L, G = 30L, force = 3,
                                 tol = 1e-10, max.iteration = 5000)
  expect_identical(out$groups, seq_along(co))              # 12 controls: one per group
  ## r.noise = 0: control j's noise is its outcome net of the refit's time
  ## effects (fitted without j), demeaned, times sqrt(T / (T - 1))
  j <- 1
  keep <- setdiff(1:N, co[j])
  f <- fect:::fect_fe(Y = Y[, keep], X = NULL, D = D[, keep], W = NULL, I = I[, keep], II = II[, keep],
                      T.on = T.on[, keep], r.cv = 0, force = 3, hasRevs = 0, tol = 1e-10,
                      max.iteration = 5000, boot = TRUE)
  R <- Y[, co[j]] - as.numeric(f$mu) - as.vector(f$xi)
  e <- (R - mean(R)) * sqrt(T / (T - 1))
  expect_equal(out$E[, j], e, tolerance = 1e-8)
})

test_that("IFE-EM bootstrap: one factor with 30 or more controls takes noise at rank 0", {
  skip_on_cran()
  f <- suppressWarnings(suppressMessages(
    fect(Y ~ D, data = simgsynth, index = c("id", "time"), method = "ife", r = 1, CV = FALSE, se = TRUE,
         vartype = "parametric", nboots = 5, seed = 1, parallel = FALSE)))
  expect_identical(f$ife.boot$r.noise, 0L)
  expect_true(all(is.finite(f$att.avg.boot)))
  expect_false(anyNA(f$ife.boot$noise))
})

test_that("IFE-EM bootstrap: settings it does not take", {
  skip_on_cran()
  base <- list(formula = Y ~ D, data = simgsynth, index = c("id", "time"), method = "ife", r = 1,
               se = TRUE, vartype = "parametric", nboots = 5, parallel = FALSE)
  expect_error(do.call(fect, c(base, list(force = "time"))), "unit")
  expect_error(do.call(fect, c(base, list(placeboTest = TRUE, placebo.period = c(-2, 0)))), "placebo")
  ## method = "fe" keeps the old gate
  expect_error(fect(Y ~ D, data = simgsynth, index = c("id", "time"), method = "fe", se = TRUE,
                    vartype = "parametric", nboots = 5, parallel = FALSE), "nevertreated")
})
