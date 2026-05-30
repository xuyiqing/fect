## ---------------------------------------------------------------
## Shared fixtures for the score-unify test family.
##
## testthat auto-loads helper-*.R files before any test-*.R, so the
## fixtures defined here (make_factor_data, ntdata, out_base) are
## visible to every test-score-*.R file. This avoids duplicating the
## ~60 lines of fixture setup across the 7 split test files.
##
## Originally lived at the top of test-score-unify.R; extracted on
## 2026-05-03 when that file was split for readability + progress
## visibility under reporter = "summary".
## ---------------------------------------------------------------

suppressWarnings(data("simdata", package = "fect"))

## DGP with factor structure and sufficient never-treated units for CV.
## N=50, TT=20, Ntr=15 => 35 never-treated units with 20 pre-treatment
## periods --- plenty for cross-validation with r up to 3.
make_factor_data <- function(N = 50, TT = 20, Ntr = 15, tau = 3.0,
                              r = 2, seed = 42) {
  set.seed(seed)
  F_mat <- matrix(rnorm(TT * r), TT, r)
  L_mat <- matrix(rnorm(N * r), N, r)
  alpha_i <- rnorm(N, 0, 1)
  xi_t <- rnorm(TT, 0, 0.5)

  T0_vec <- rep(Inf, N)
  if (Ntr > 0) {
    T0_vec[1:Ntr] <- sample(round(TT * 0.4):round(TT * 0.7), Ntr,
                             replace = TRUE)
  }

  Y_vec <- D_vec <- numeric(N * TT)
  id_vec <- time_vec <- integer(N * TT)
  idx <- 1
  for (i in 1:N) {
    for (t in 1:TT) {
      treated <- (t >= T0_vec[i])
      D_vec[idx] <- as.integer(treated)
      Y_vec[idx] <- alpha_i[i] + xi_t[t] +
        sum(F_mat[t, ] * L_mat[i, ]) +
        tau * D_vec[idx] + rnorm(1, 0, 0.5)
      id_vec[idx] <- i
      time_vec[idx] <- t
      idx <- idx + 1
    }
  }

  data.frame(id = id_vec, time = time_vec, Y = Y_vec, D = D_vec)
}

## Shared never-treated fixture for Section C / E / F tests.
ntdata <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

## Shared fitted object for the .score_residuals() / fect_cv / fect_mspe
## test families. Fit once at session start; all test-*.R files that need
## it can reference the global symbol.
out_base <- suppressWarnings(suppressMessages(
  fect::fect(
    Y ~ D + X1 + X2,
    data    = simdata,
    index   = c("id", "time"),
    method  = "ife",
    r       = 2,
    CV      = FALSE,
    se      = FALSE,
    parallel = FALSE
  )
))
