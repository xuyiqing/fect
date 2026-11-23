## The exact inversion of full conformal inference finds where two lines cross
## by dividing by the difference of their slopes. The slopes come from a
## difference quotient, so two slopes that are equal (or opposite) in exact
## arithmetic can differ in their last digits. Before 2.5.1 that gave a root
## near 1e14 times the errors, a set with an extra piece beyond it, and a hull
## that lost a finite end. A factor's slope that is zero up to rounding is
## now taken to be zero (.fc_factors in R/conformal-full.R).

fcrr_coefs <- function() {
  ## unit 1 is the treated unit (slope -1); unit 2's prediction puts weight 1
  ## on the treated unit (slope 1); the others put little or none
  list(a = c(0.5, 0.3, -1.2, 2.0, 0.7, -0.4, 1.5, -2.2, 0.9, 1.1),
       b = c(-1, 1, 0, 0, 0.2, 0, 0, 0, 0.1, 0))
}

## With 10 units and gamma = 0.15 a candidate is in the set as soon as one
## control's score is at least the treated unit's, so one comparison that
## flips far out changes the set: with a slope just above 1 the set gains a
## piece on the far left, with a slope just below 1 it is cut off on the far
## right.
test_that("a slope that differs from 1 in its last digit gives the same set as 1", {
  cf <- fcrr_coefs()
  exact <- fect:::.fc_invert_affine(cf$a, cf$b, 1L, 0.15)
  for (eps in c(1, -1, 4, -4) * .Machine$double.eps) {
    b <- cf$b
    b[2] <- 1 + eps
    expect_false(b[2] == 1)
    expect_identical(fect:::.fc_invert_affine(cf$a, b, 1L, 0.15), exact)
  }
  ## no end of the set is a rounding root, and the set is unbounded above
  ends <- c(exact$components$lower, exact$components$upper)
  expect_true(all(abs(ends[is.finite(ends)]) < 1e3))
  expect_identical(exact$CI.upper, Inf)
  expect_true(is.finite(exact$CI.lower))
})

test_that("the same holds for the studentized score, where both slopes are rescaled", {
  cf <- fcrr_coefs()
  w <- c(2e-6, 2e-6, 3e-6, 1e-6, 4e-6, 2e-6, 5e-6, 1e-6, 3e-6, 2e-6)
  exact <- fect:::.fc_invert_affine(cf$a * w, cf$b * w, 1L, 0.15)
  b <- cf$b
  b[2] <- 1 + .Machine$double.eps
  expect_identical(fect:::.fc_invert_affine(cf$a * w, b * w, 1L, 0.15), exact)
  ## a slope sum that is small but not rounding (1e-6 of the slopes) keeps its root
  b[2] <- 1 + 1e-6
  far <- fect:::.fc_root_table(cf$a * w, b * w, 1L)
  expect_true(any(abs(far$root) > 1e5))
})

test_that("a set that is a union of intervals in exact arithmetic stays one", {
  ## one control with slope 2: |2c - 10| >= |c| holds for c <= 10/3 and c >= 10
  v <- fect:::.fc_invert_affine(c(0, -10), c(-1, 2), 1L, 0.6)
  expect_true(v$disconnected)
  expect_equal(nrow(v$components), 2L)
  expect_equal(v$components$upper[1], 10 / 3)
  expect_equal(v$components$lower[2], 10)
  expect_identical(c(v$CI.lower, v$CI.upper), c(-Inf, Inf))
})

test_that("conformal_sc() with a near-twin control reports a one-sided set, not the whole line", {
  skip_on_cran()
  skip_if_not_installed("quadprog")
  ## one treated unit (u01), nine controls, 24 dates, adoption at date 19; two
  ## controls are near copies of the treated unit, so their synthetic controls
  ## put all, or nearly all, of their weight on it
  set.seed(445)
  N <- 10; TT <- 24
  alpha <- rnorm(N); xi <- rnorm(TT)
  lam <- matrix(rnorm(N * 2), N, 2); f <- matrix(rnorm(TT * 2), TT, 2)
  lam[1, ] <- lam[1, ] + 0.5
  x1 <- matrix(rnorm(N * TT), TT, N)
  eps <- matrix(rnorm(N * TT, sd = 1), TT, N)
  Y <- matrix(xi, TT, N) + matrix(alpha, TT, N, byrow = TRUE) + f %*% t(lam) + x1 + eps
  D <- matrix(0L, TT, N); D[19:TT, 1] <- 1L
  Y <- Y + 2 * D
  set.seed(5445)
  y1 <- Y[, 1] - 2 * D[, 1]
  for (j in 2:3) Y[, j] <- y1 + rnorm(TT, sd = 0.3)
  d <- data.frame(id = rep(sprintf("u%02d", 1:N), each = TT), time = rep(2000L + 1:TT, N),
                  Y = as.vector(Y), D = as.vector(D))
  msg <- character(0)
  fit <- withCallingHandlers(
    conformal_sc(Y ~ D, data = d, index = c("id", "time"), rule = "sc", alpha = 0.1, diag = 4),
    message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") })
  rows <- rbind(fit$unit, fit$diag)
  expect_false(any(rows$disconnected %in% TRUE))
  expect_false(any(grepl("unions of intervals", msg)))
  ## no end of any set is a rounding root
  ends <- unlist(lapply(c(fit$components, fit$diag.components),
                        function(cp) c(cp$lower, cp$upper)))
  expect_true(all(abs(ends[is.finite(ends)]) < 1e6))
  ## the date whose set was reported as the whole line has a finite upper end
  r20 <- fit$unit[fit$unit$time == 2020, ]
  expect_identical(r20$CI.lower, -Inf)
  expect_true(is.finite(r20$CI.upper))
})
