## Validity of full conformal inference by simulation (Liu and Xu 2026,
## coverage theorem): with exchangeable units the set covers the effect with
## probability at least 1 - floor(gamma N)/N, and the p-value at the true
## effect is super-uniform on {1/N, ..., 1}.

fc_exch_draw <- function(J, T0, tau) {
  N <- J + 1
  TT <- T0 + 1
  alpha.i <- rnorm(N)
  xi <- cumsum(rnorm(TT, 0, 0.5))
  E <- matrix(rnorm(TT * N), TT, N)
  Y <- outer(xi, rep(1, N)) + matrix(alpha.i, TT, N, byrow = TRUE) + E
  Y[TT, 1] <- Y[TT, 1] + tau
  list(Ypre = Y[1:T0, , drop = FALSE], ytg = Y[TT, ])
}

test_that("full conformal sets cover at the guaranteed rate (DID, exchangeable units)", {
  skip_on_cran()
  set.seed(20261111)
  J <- 19
  gamma <- 0.1
  tau <- 1.5
  R <- 2000
  cover <- logical(R)
  pv <- numeric(R)
  for (r in seq_len(R)) {
    dd <- fc_exch_draw(J, T0 = 8, tau = tau)
    cf <- fect:::.fc_coefs(dd$Ypre, dd$ytg, 1L, fect:::.fc_rule_did)
    v <- fect:::.fc_invert_affine(cf$a, cf$b, 1L, gamma)
    comp <- v$components
    cover[r] <- any(comp$lower <= tau & tau <= comp$upper)
    pv[r] <- fect:::.fc_affine_p(cf$a, cf$b, 1L, tau)
  }
  N <- J + 1
  level <- 1 - floor(gamma * N) / N            # 0.9
  se <- sqrt(level * (1 - level) / R)
  expect_gte(mean(cover), level - 3 * se)
  ## super-uniform p-values at the truth: P(p <= k/N) <= k/N
  for (k in 1:5) {
    expect_lte(mean(pv <= k / N), k / N + 3 * sqrt((k / N) * (1 - k / N) / R))
  }
})

test_that("the Jackknife+ covers at its guaranteed rate (DID, exchangeable units)", {
  skip_on_cran()
  set.seed(20261112)
  J <- 19
  gamma <- 0.1
  tau <- -1
  R <- 1500
  cover <- logical(R)
  for (r in seq_len(R)) {
    dd <- fc_exch_draw(J, T0 = 8, tau = tau)
    v <- fect:::.fc_jplus(dd$Ypre, dd$ytg, 1L, fect:::.fc_rule_did, gamma)
    cover[r] <- v$CI.lower <= tau && tau <= v$CI.upper
  }
  N <- J + 1
  q <- fect:::.conformal_q(N, gamma)
  level <- 1 - (2 * q - 1) / N                 # 0.85 at N = 20
  expect_gte(mean(cover), level - 3 * sqrt(level * (1 - level) / R))
})
