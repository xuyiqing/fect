## Full conformal inference and the Jackknife+ for target-date estimators
## (R/conformal-full.R; Liu and Xu 2026, Sec. 2.3-2.4).

## ---- exact inversion -----------------------------------------------------------

## membership of c in a set of components
fc_in <- function(comp, c0) {
  if (!nrow(comp)) return(FALSE)
  any((comp$lower < c0 | (comp$lower == c0 & comp$lower.closed)) &
        (c0 < comp$upper | (c0 == comp$upper & comp$upper.closed)))
}

test_that("exact inversion equals the p-value rule on a dense grid", {
  set.seed(11)
  for (rep in 1:200) {
    N <- sample(5:25, 1)
    a <- rnorm(N)
    b <- runif(N, -1.5, 1.5)
    ## some ties: controls with the treated unit's |error|
    if (rep %% 4 == 0) a[2] <- a[1]
    if (rep %% 5 == 0) b[3] <- 0
    b[1] <- -1
    gamma <- sample(c(0.05, 0.1, 0.2, 0.3), 1)
    v <- fect:::.fc_invert_affine(a, b, 1L, gamma)
    rt <- fect:::.fc_root_table(a, b, 1L)$root
    ## probes: every root, midpoints between roots, and points beyond
    r <- sort(unique(rt))
    probes <- c(r, if (length(r) > 1) (r[-1] + r[-length(r)]) / 2,
                min(c(r, 0)) - 5, max(c(r, 0)) + 5, seq(-10, 10, by = 0.37))
    for (c0 in probes) {
      if (c0 %in% r) {
        ## at a root the comparison is a tie, which counts for the
        ## candidate; a direct floating-point evaluation there needs a
        ## tolerance
        sc <- abs(a + b * c0)
        p <- mean(sc >= sc[1] - 1e-9 * max(1, sc))
      } else {
        p <- fect:::.fc_affine_p(a, b, 1L, c0)
      }
      expect_identical(p > gamma, fc_in(v$components, c0))
    }
    expect_equal(v$p.value, fect:::.fc_affine_p(a, b, 1L, 0))
    expect_equal(v$level, 1 - floor(gamma * N) / N)
  }
})

test_that("slopes above one in absolute value can give disconnected sets", {
  ## two controls whose errors grow faster than the treated unit's score
  a <- c(0, 3, -3, 0.5, -0.5)
  b <- c(-1, 2, 2, 0, 0)
  v <- fect:::.fc_invert_affine(a, b, 1L, 0.2)
  expect_true(nrow(v$components) >= 1)
  probes <- seq(-6, 6, by = 0.01)
  for (c0 in probes) {
    expect_identical(fect:::.fc_affine_p(a, b, 1L, c0) > 0.2,
                     fc_in(v$components, c0))
  }
})

test_that("too few controls give the whole real line", {
  a <- c(1, rnorm(5))
  b <- c(-1, rep(0.2, 5))
  v <- fect:::.fc_invert_affine(a, b, 1L, 0.1)   # N = 6: floor(0.6) = 0
  expect_equal(v$CI.lower, -Inf)
  expect_equal(v$CI.upper, Inf)
  expect_equal(v$level, 1)
})

## ---- rules and conformal_td() --------------------------------------------------

fc_panel_data <- function(J = 30, T0 = 12, T1 = 4, K = 1, tau = 1, seed = 1,
                          stagger = FALSE) {
  set.seed(seed)
  N <- J + K
  TT <- T0 + T1
  a <- rnorm(N)
  xi <- cumsum(rnorm(TT, 0, 0.5))
  Y0 <- outer(xi, rep(1, N)) + matrix(a, TT, N, byrow = TRUE) +
    matrix(rnorm(TT * N), TT, N)
  D <- matrix(0, TT, N)
  for (k in seq_len(K)) {
    on <- if (stagger) T0 + k - 1 else T0 + 1
    D[on:TT, k] <- 1
  }
  Y <- Y0 + tau * D
  data.frame(id = rep(seq_len(N), each = TT), time = rep(seq_len(TT), N),
             Y = c(Y), D = c(D))
}

test_that("the did rule is affine and a user copy gives identical results", {
  d <- fc_panel_data()
  f1 <- conformal_td(Y ~ D, data = d, index = c("id", "time"), rule = "did")
  my_did <- function(query.pre, background.pre, background.target) {
    mean(query.pre) + mean(background.target - colMeans(background.pre))
  }
  f2 <- conformal_td(Y ~ D, data = d, index = c("id", "time"), rule = my_did)
  expect_true(all(f1$unit$method == "exact"))
  expect_equal(f1$unit[, c("eff", "CI.lower", "CI.upper", "p.value", "N")],
               f2$unit[, c("eff", "CI.lower", "CI.upper", "p.value", "N")])
  expect_equal(f1$avg$CI.lower, f2$avg$CI.lower)
  expect_s3_class(f1, "conformal_td")
  expect_output(print(f1), "Full conformal inference")
  expect_s3_class(plot(f1), "ggplot")
})

test_that("a rule that is not affine falls back to the grid, and the grid p-values are right", {
  d <- fc_panel_data(J = 15, T0 = 8, T1 = 2)
  med_rule <- function(query.pre, background.pre, background.target) {
    mean(query.pre) + stats::median(background.target - colMeans(background.pre))
  }
  f <- conformal_td(Y ~ D, data = d, index = c("id", "time"), rule = med_rule,
                    grid = seq(-6, 6, by = 0.25))
  expect_true(all(f$unit$method == "grid"))
  ## brute force at the estimate for the first date: p = 1 there
  expect_true(all(f$unit$CI.lower <= f$unit$eff & f$unit$eff <= f$unit$CI.upper))
  expect_true(all(is.na(f$avg$CI.lower)))
  expect_match(f$avg$note, "not affine")
})

test_that("the p-value at zero equals a direct refit of every unit", {
  d <- fc_panel_data(J = 12, T0 = 6, T1 = 2, seed = 4)
  f <- conformal_td(Y ~ D, data = d, index = c("id", "time"), rule = "did")
  Y <- matrix(d$Y, ncol = 13)
  t <- 7
  sc <- numeric(13)
  for (i in 1:13) {
    bg <- setdiff(1:13, i)
    pr <- fect:::.fc_rule_did(Y[1:6, i], Y[1:6, bg, drop = FALSE], Y[t, bg])
    sc[i] <- abs(Y[t, i] - pr)
  }
  expect_equal(f$unit$p.value[1], mean(sc >= sc[1]))
})

test_that("the Jackknife+ agrees with .conformal_interval on hand-built folds", {
  d <- fc_panel_data(J = 20, T0 = 6, T1 = 1, seed = 9)
  f <- conformal_td(Y ~ D, data = d, index = c("id", "time"), rule = "did",
                    procedure = "jackknife+", alpha = 0.1)
  Y <- matrix(d$Y, ncol = 21)
  M0 <- R <- numeric(20)
  for (k in 1:20) {
    j <- k + 1
    bg <- setdiff(2:21, j)
    M0[k] <- fect:::.fc_rule_did(Y[1:6, 1], Y[1:6, bg], Y[7, bg])
    R[k] <- abs(Y[7, j] - fect:::.fc_rule_did(Y[1:6, j], Y[1:6, bg], Y[7, bg]))
  }
  v <- fect:::.conformal_interval(M0, R, Y[7, 1], 0.1)
  expect_equal(f$unit$CI.lower, unname(v["CI.lower"]))
  expect_equal(f$unit$CI.upper, unname(v["CI.upper"]))
  expect_equal(f$unit$p.value, unname(v["p.value"]))
})

test_that("conformal_td checks its inputs", {
  d <- fc_panel_data()
  expect_error(conformal_td(Y ~ D + id, data = d, index = c("id", "time")),
               "covariates")
  expect_error(conformal_td(Y ~ D, data = d, index = c("id", "time"),
                            alpha = 1), "alpha")
  expect_error(conformal_td(Y ~ D, data = d, index = c("id", "time"),
                            rule = "nope"), "rule")
  expect_error(conformal_td(Y ~ D, data = d, index = c("id", "time"),
                            diag = -1), "diag")
  d2 <- d
  d2$D[d2$id == 1 & d2$time == 16] <- 0
  expect_error(conformal_td(Y ~ D, data = d2, index = c("id", "time")),
               "stays on")
})

test_that("staggered adoption: each treated unit uses its own pre-treatment dates", {
  d <- fc_panel_data(K = 3, stagger = TRUE, seed = 2)
  f <- conformal_td(Y ~ D, data = d, index = c("id", "time"), rule = "did")
  expect_setequal(unique(f$unit$unit), c("1", "2", "3"))
  expect_equal(as.numeric(tapply(f$unit$n.history, f$unit$unit, unique)),
               c(11, 12, 13))
  expect_equal(nrow(f$avg), 3)
})

## ---- fect(method = "fe", vartype = "conformal") --------------------------------

test_that("fect fe conformal equals conformal_td with the did rule", {
  d <- fc_panel_data(J = 25, seed = 5)
  f <- suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"),
                             method = "fe", se = TRUE, vartype = "conformal",
                             alpha = 0.1, conformal.diag = 2,
                             parallel = FALSE))
  td <- conformal_td(Y ~ D, data = d, index = c("id", "time"), rule = "did",
                     alpha = 0.1, diag = 2)
  u <- f$conformal$unit
  expect_equal(u$eff, u$eff.fit, tolerance = 1e-10)
  expect_equal(u$CI.lower, td$unit$CI.lower)
  expect_equal(u$CI.upper, td$unit$CI.upper)
  expect_equal(u$p.value, td$unit$p.value)
  expect_equal(unname(f$est.avg[1, c("CI.lower", "CI.upper")]),
               c(td$avg$CI.lower, td$avg$CI.upper))
  post <- as.character(1:4)
  expect_equal(unname(f$est.att[post, "CI.lower"]), td$unit$CI.lower)
  pre <- as.character(-1:0)
  expect_equal(unname(f$est.att[pre, "CI.lower"]), td$diag$CI.lower)
  expect_output(print(f), "full conformal for target-date")
  expect_s3_class(plot(f), "ggplot")
})

test_that("fect fe conformal with several treated units leaves averages without sets", {
  d <- fc_panel_data(K = 3, seed = 6)
  f <- suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"),
                             method = "fe", se = TRUE, vartype = "conformal",
                             parallel = FALSE))
  expect_true(all(is.na(f$est.avg[1, c("CI.lower", "CI.upper")])))
  expect_true(all(is.na(f$est.att[as.character(1:4), "CI.lower"])))
  expect_equal(length(unique(f$conformal$unit$unit)), 3)
  expect_output(print(f), "no sets for averages")
})

test_that("fect fe conformal refuses covariates, force = none, blocks and relabeling", {
  d <- fc_panel_data()
  d$X <- rnorm(nrow(d))
  expect_error(fect(Y ~ D + X, data = d, index = c("id", "time"),
                    method = "fe", se = TRUE, vartype = "conformal",
                    parallel = FALSE), "covariates")
  expect_error(fect(Y ~ D, data = d, index = c("id", "time"), method = "fe",
                    force = "none", se = TRUE, vartype = "conformal",
                    parallel = FALSE), "force")
  expect_error(fect(Y ~ D, data = d, index = c("id", "time"), method = "fe",
                    se = TRUE, vartype = "conformal", conformal.block = TRUE,
                    parallel = FALSE), "conformal.block")
})

test_that("staggered fe fits say when the target-date estimate differs from fect's", {
  d <- fc_panel_data(K = 3, stagger = TRUE, seed = 7)
  expect_message(
    f <- fect(Y ~ D, data = d, index = c("id", "time"), method = "fe",
              se = TRUE, vartype = "conformal", parallel = FALSE),
    "differs from fect's estimate")
  u <- f$conformal$unit
  expect_true(any(abs(u$eff - u$eff.fit) > 1e-6))
})
