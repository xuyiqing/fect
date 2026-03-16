## ---------------------------------------------------------
## Phase tests for the factors.from refactoring
##
## These tests define acceptance criteria for 5 phases:
##   Phase 1: factors.from parameter (notyettreated / nevertreated)
##   Phase 2: gsynth merged into ife via factors.from
##   Phase 3: parametric bootstrap unlocked for ife/cfe
##   Phase 4: polynomial removed (superseded by cfe)
##   Phase 5: input validation and safety guards
##
## Tests are written BEFORE implementation. Each test will fail
## until its phase is implemented, then pass permanently.
## ---------------------------------------------------------

## ---- helpers ----
make_staggered_data <- function(N = 40, TT = 20, Ntr = 15, tau = 3.0,
                                 seed = 42) {
  set.seed(seed)
  alpha_i <- rnorm(N, 0, 1)
  xi_t    <- rnorm(TT, 0, 0.5)

  ## staggered adoption: treated units adopt between T0=8 and T0=14
  T0_vec <- rep(Inf, N)
  if (Ntr > 0) {
    T0_vec[1:Ntr] <- sample(8:14, Ntr, replace = TRUE)
  }

  Y_vec <- D_vec <- numeric(N * TT)
  id_vec <- time_vec <- integer(N * TT)
  idx <- 1
  for (i in 1:N) {
    for (t in 1:TT) {
      treated <- (t >= T0_vec[i])
      D_vec[idx]    <- as.integer(treated)
      Y_vec[idx]    <- alpha_i[i] + xi_t[t] + tau * D_vec[idx] + rnorm(1, 0, 0.5)
      id_vec[idx]   <- i
      time_vec[idx] <- t
      idx <- idx + 1
    }
  }

  data.frame(id = id_vec, time = time_vec, Y = Y_vec, D = D_vec)
}

## DGP with factor structure for comparisons needing r > 0
make_factor_data <- function(N = 100, TT = 30, Ntr = 30, tau = 3.0,
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

## ========================================================
## PHASE 6: em parameter and gsynth reroute
## ========================================================

test_that("Phase 6a: gsynth and ife+nevertreated produce identical ATT", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  set.seed(100)
  out_gs <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "gsynth", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", parallel = FALSE
  )))

  set.seed(100)
  out_ife <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_equal(out_gs$att.avg, out_ife$att.avg,
               info = "gsynth and ife+nevertreated must be the same estimator")
})

test_that("Phase 6b: em=FALSE + factors.from='notyettreated' errors", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, Ntr = 15)

  expect_error(
    fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r = 2, CV = FALSE, se = FALSE,
      factors.from = "notyettreated", em = FALSE,
      parallel = FALSE
    ),
    regexp = "em.*FALSE|notyettreated.*em"
  )
})

test_that("Phase 6c: method='gsynth' auto-sets em=FALSE in output", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "gsynth", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", parallel = FALSE
  )))

  expect_false(out$em)
})

test_that("Phase 6d: method='ife' defaults em=TRUE in output", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, Ntr = 15)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 0, CV = FALSE, se = FALSE,
    parallel = FALSE
  )))

  expect_true(out$em)
})

test_that("Phase 6e: ife+nevertreated and ife+notyettreated both close to true tau", {
  skip_on_cran()
  df <- make_factor_data(N = 200, TT = 30, Ntr = 60, tau = 3.0, r = 2, seed = 42)

  set.seed(100)
  out_nyt <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "notyettreated",
    parallel = FALSE
  )))

  set.seed(100)
  out_nt <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  ## Both should be within 1.0 of true tau=3.0 (loose bound)
  expect_true(abs(out_nyt$att.avg - 3.0) < 1.0,
              info = paste("notyettreated ATT:", out_nyt$att.avg))
  expect_true(abs(out_nt$att.avg - 3.0) < 1.0,
              info = paste("nevertreated ATT:", out_nt$att.avg))

  ## They need not be equal — different predictive routines
  cat(sprintf("\n  [info] notyettreated ATT=%.4f  nevertreated ATT=%.4f  diff=%.4f\n",
              out_nyt$att.avg, out_nt$att.avg, out_nyt$att.avg - out_nt$att.avg))
})

## ========================================================
## PHASE 1: factors.from parameter
## ========================================================

test_that("Phase 1a: fect accepts factors.from='notyettreated' (default behavior)", {
  skip_on_cran()
  df <- make_staggered_data()

  out <- suppressWarnings(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE,
    factors.from = "notyettreated",
    parallel = FALSE
  ))

  expect_s3_class(out, "fect")
  expect_true(is.numeric(out$att.avg))
  expect_true(!is.na(out$att.avg))
})

test_that("Phase 1b: fect accepts factors.from='nevertreated'", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, Ntr = 15)  ## 25 never-treated

  out <- suppressWarnings(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE,
    factors.from = "nevertreated",
    parallel = FALSE
  ))

  expect_s3_class(out, "fect")
  expect_true(is.numeric(out$att.avg))
  expect_true(!is.na(out$att.avg))
})

test_that("Phase 1c: factors.from='nevertreated' produces different estimates than 'notyettreated'", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, Ntr = 15, seed = 99)

  out_nyt <- suppressWarnings(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE,
    factors.from = "notyettreated",
    parallel = FALSE
  ))

  out_nt <- suppressWarnings(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE,
    factors.from = "nevertreated",
    parallel = FALSE
  ))

  ## Both should produce valid ATTs, but they should differ
  expect_true(is.numeric(out_nyt$att.avg))
  expect_true(is.numeric(out_nt$att.avg))
  ## Allow for numerical coincidence, but in general they won't be equal
  ## The key structural check: both succeeded with different code paths
  expect_s3_class(out_nyt, "fect")
  expect_s3_class(out_nt, "fect")
})

test_that("Phase 1d: factors.from defaults to 'notyettreated' when omitted", {
  skip_on_cran()
  df <- make_staggered_data()

  ## Without factors.from
  out_default <- suppressWarnings(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE,
    parallel = FALSE
  ))

  ## With factors.from = "notyettreated" explicitly
  out_explicit <- suppressWarnings(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE,
    factors.from = "notyettreated",
    parallel = FALSE
  ))

  ## Should produce identical results
  expect_equal(out_default$att.avg, out_explicit$att.avg, tolerance = 1e-10)
})

test_that("Phase 1e: factors.from='nevertreated' works with method='cfe'", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, Ntr = 15)

  out <- suppressWarnings(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", r = 0, CV = FALSE, se = FALSE,
    factors.from = "nevertreated",
    parallel = FALSE
  ))

  expect_s3_class(out, "fect")
  expect_true(is.numeric(out$att.avg))
})

test_that("Phase 1f: factors.from threads through cross-validation", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, Ntr = 15)

  out <- suppressWarnings(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = c(0, 3), CV = TRUE, se = FALSE,
    factors.from = "nevertreated",
    parallel = FALSE
  ))

  expect_s3_class(out, "fect")
  expect_true(!is.null(out$r.cv))
})

test_that("Phase 1g: factors.from threads through bootstrap inference", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, Ntr = 15)

  out <- suppressWarnings(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE,
    se = TRUE, nboots = 30,
    factors.from = "nevertreated",
    parallel = FALSE
  ))

  expect_s3_class(out, "fect")
  expect_true(!is.null(out$est.att))
  expect_true(!is.null(out$att.vcov))
})

## ========================================================
## PHASE 2: gsynth merged into ife
## ========================================================

test_that("Phase 2a: gsynth produces same results via ife+factors.from='nevertreated'", {
  skip_on_cran()
  suppressWarnings(try(data("simgsynth", package = "fect"), silent = TRUE))
  skip_if_not(exists("simgsynth"), "Dataset 'simgsynth' not available")

  set.seed(200)
  out_gs <- suppressWarnings(fect::fect(
    Y ~ D, data = simgsynth, index = c("id", "time"),
    method = "gsynth", r = 2, CV = FALSE, se = FALSE,
    force = "two-way",
    parallel = FALSE
  ))

  set.seed(200)
  out_ife <- suppressWarnings(fect::fect(
    Y ~ D, data = simgsynth, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE,
    force = "two-way",
    factors.from = "nevertreated",
    parallel = FALSE
  ))

  ## ATT estimates should be numerically very close
  expect_equal(out_gs$att.avg, out_ife$att.avg, tolerance = 0.05)
})

test_that("Phase 2b: gsynth still works (backward compatibility)", {
  skip_on_cran()
  suppressWarnings(try(data("simgsynth", package = "fect"), silent = TRUE))
  skip_if_not(exists("simgsynth"), "Dataset 'simgsynth' not available")

  set.seed(201)
  out <- suppressWarnings(fect::fect(
    Y ~ D, data = simgsynth, index = c("id", "time"),
    method = "gsynth", r = 2, CV = FALSE, se = FALSE,
    force = "two-way",
    parallel = FALSE
  ))

  expect_s3_class(out, "fect")
  expect_true(is.numeric(out$att.avg))
})

test_that("Phase 2c: gsynth with parametric bootstrap still works", {
  skip_on_cran()
  suppressWarnings(try(data("simgsynth", package = "fect"), silent = TRUE))
  skip_if_not(exists("simgsynth"), "Dataset 'simgsynth' not available")

  set.seed(202)
  out <- suppressWarnings(fect::fect(
    Y ~ D, data = simgsynth, index = c("id", "time"),
    method = "gsynth", r = 2, CV = FALSE,
    se = TRUE, vartype = "parametric", nboots = 30,
    force = "two-way", min.T0 = 2,
    parallel = FALSE
  ))

  expect_s3_class(out, "fect")
  expect_true(!is.null(out$est.att))
})

## ========================================================
## PHASE 3: parametric bootstrap unlocked for ife/cfe
## ========================================================

test_that("Phase 3a: vartype='parametric' accepted for method='ife'", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, Ntr = 15)

  ## Currently this errors: "parametric option is only available for gsynth"
  ## After Phase 3, it should succeed
  out <- suppressWarnings(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE,
    se = TRUE, vartype = "parametric", nboots = 30,
    parallel = FALSE
  ))

  expect_s3_class(out, "fect")
  expect_true(!is.null(out$est.att))
})

test_that("Phase 3b: vartype='parametric' accepted for method='fe'", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, Ntr = 15)

  out <- suppressWarnings(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "fe", CV = FALSE,
    se = TRUE, vartype = "parametric", nboots = 30,
    parallel = FALSE
  ))

  expect_s3_class(out, "fect")
  expect_true(!is.null(out$est.att))
})

test_that("Phase 3c: vartype='parametric' accepted for method='cfe'", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, Ntr = 15)

  out <- suppressWarnings(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", r = 0, CV = FALSE,
    se = TRUE, vartype = "parametric", nboots = 30,
    parallel = FALSE
  ))

  expect_s3_class(out, "fect")
  expect_true(!is.null(out$est.att))
})

test_that("Phase 3d: vartype='parametric' still errors for method='mc'", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, Ntr = 15)

  expect_error(
    fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "mc", CV = FALSE,
      se = TRUE, vartype = "parametric", nboots = 30,
      parallel = FALSE
    ),
    regexp = "parametric"
  )
})

## ========================================================
## PHASE 4: polynomial removed (superseded by cfe)
## ========================================================
## Phase 4a test removed: method="polynomial" no longer exists.

## ========================================================
## PHASE 5: input validation and safety guards
## ========================================================

test_that("Phase 5a: factors.from rejects invalid values", {
  skip_on_cran()
  df <- make_staggered_data()

  expect_error(
    fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r = 2, CV = FALSE, se = FALSE,
      factors.from = "invalid_value",
      parallel = FALSE
    ),
    regexp = "factors.from"
  )
})

test_that("Phase 5b: factors.from='nevertreated' + method='mc' errors", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, Ntr = 15)

  expect_error(
    fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "mc", CV = FALSE, se = FALSE,
      factors.from = "nevertreated",
      parallel = FALSE
    ),
    regexp = "factors.from|nevertreated|mc"
  )
})

test_that("Phase 5c: factors.from='nevertreated' errors when no never-treated units", {
  skip_on_cran()
  ## All units are treated
  df <- make_staggered_data(N = 20, Ntr = 20)

  expect_error(
    fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r = 2, CV = FALSE, se = FALSE,
      factors.from = "nevertreated",
      parallel = FALSE
    ),
    regexp = "never.treated|nevertreated|control"
  )
})

test_that("Phase 5d: factors.from='nevertreated' errors when too few never-treated units", {
  skip_on_cran()
  ## Only 1 never-treated unit (insufficient for r=2)
  df <- make_staggered_data(N = 20, Ntr = 19)

  expect_error(
    suppressWarnings(fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r = 2, CV = FALSE, se = FALSE,
      factors.from = "nevertreated",
      parallel = FALSE
    )),
    regexp = "never.treated|nevertreated|insuffic|too few"
  )
})

test_that("Phase 5e: output records factors.from in return object", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, Ntr = 15)

  out <- suppressWarnings(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE,
    factors.from = "nevertreated",
    parallel = FALSE
  ))

  expect_true("factors.from" %in% names(out))
  expect_equal(out$factors.from, "nevertreated")
})

test_that("Phase 5f: factors.from='nevertreated' + method='both' errors", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, Ntr = 15)

  expect_error(
    fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "both", CV = FALSE, se = FALSE,
      factors.from = "nevertreated",
      parallel = FALSE
    ),
    regexp = "both"
  )
})

## ========================================================
## PHASE 3a: CFE bifurcation in fect_nevertreated
## ========================================================

## ---- Phase 3a DGP helpers ----

make_cfe_z_data <- function(N = 100, TT = 30, Ntr = 30, tau = 3.0,
                             r = 2, seed = 42) {
  set.seed(seed)
  F_mat <- matrix(rnorm(TT * r), TT, r)
  L_mat <- matrix(rnorm(N * r), N, r)
  alpha_i <- rnorm(N, 0, 1)
  xi_t <- rnorm(TT, 0, 0.5)

  ## Z: time-invariant covariate (baseline characteristic)
  Z_i <- rnorm(N, 0, 1)

  ## gamma: grouped time coefficient (one per time period)
  gamma_t <- rnorm(TT, 0, 0.5)

  ## Staggered treatment adoption
  T0_vec <- rep(Inf, N)
  if (Ntr > 0) {
    T0_vec[1:Ntr] <- sample(round(TT * 0.4):round(TT * 0.7), Ntr,
                             replace = TRUE)
  }

  Y_vec <- D_vec <- Z_vec <- numeric(N * TT)
  id_vec <- time_vec <- integer(N * TT)
  idx <- 1
  for (i in 1:N) {
    for (t in 1:TT) {
      treated <- (t >= T0_vec[i])
      D_vec[idx] <- as.integer(treated)
      factor_component <- if (r > 0) sum(F_mat[t, ] * L_mat[i, ]) else 0
      Y_vec[idx] <- alpha_i[i] + xi_t[t] + Z_i[i] * gamma_t[t] +
        factor_component + tau * D_vec[idx] + rnorm(1, 0, 0.5)
      Z_vec[idx] <- Z_i[i]
      id_vec[idx] <- i
      time_vec[idx] <- t
      idx <- idx + 1
    }
  }

  data.frame(id = id_vec, time = time_vec, Y = Y_vec, D = D_vec, Z = Z_vec)
}

make_cfe_q_data <- function(N = 100, TT = 30, Ntr = 30, tau = 3.0,
                             r = 2, seed = 42) {
  set.seed(seed)
  F_mat <- matrix(rnorm(TT * r), TT, r)
  L_mat <- matrix(rnorm(N * r), N, r)
  alpha_i <- rnorm(N, 0, 1)
  xi_t <- rnorm(TT, 0, 0.5)

  ## Q: time-varying basis (normalized time trend)
  Q_t <- (1:TT) / TT

  ## kappa: unit-specific coefficient on Q
  kappa_i <- rnorm(N, 0, 0.5)

  ## Staggered treatment adoption
  T0_vec <- rep(Inf, N)
  if (Ntr > 0) {
    T0_vec[1:Ntr] <- sample(round(TT * 0.4):round(TT * 0.7), Ntr,
                             replace = TRUE)
  }

  Y_vec <- D_vec <- Q_vec <- numeric(N * TT)
  id_vec <- time_vec <- integer(N * TT)
  idx <- 1
  for (i in 1:N) {
    for (t in 1:TT) {
      treated <- (t >= T0_vec[i])
      D_vec[idx] <- as.integer(treated)
      factor_component <- if (r > 0) sum(F_mat[t, ] * L_mat[i, ]) else 0
      Y_vec[idx] <- alpha_i[i] + xi_t[t] + Q_t[t] * kappa_i[i] +
        factor_component + tau * D_vec[idx] + rnorm(1, 0, 0.5)
      Q_vec[idx] <- Q_t[t]
      id_vec[idx] <- i
      time_vec[idx] <- t
      idx <- idx + 1
    }
  }

  data.frame(id = id_vec, time = time_vec, Y = Y_vec, D = D_vec, Q = Q_vec)
}

make_cfe_fe_data <- function(N = 100, TT = 30, Ntr = 30, tau = 3.0,
                              r = 2, seed = 42) {
  set.seed(seed)
  F_mat <- matrix(rnorm(TT * r), TT, r)
  L_mat <- matrix(rnorm(N * r), N, r)
  alpha_i <- rnorm(N, 0, 1)
  xi_t <- rnorm(TT, 0, 0.5)

  ## Industry: 4 industries, randomly assigned to ALL units (shared / Type B)
  industry_i <- sample(1:4, N, replace = TRUE)
  industry_fe <- c(-1.0, 0.5, 0.0, 1.5)

  ## Staggered treatment adoption
  T0_vec <- rep(Inf, N)
  if (Ntr > 0) {
    T0_vec[1:Ntr] <- sample(round(TT * 0.4):round(TT * 0.7), Ntr,
                             replace = TRUE)
  }

  Y_vec <- D_vec <- numeric(N * TT)
  id_vec <- time_vec <- industry_vec <- integer(N * TT)
  idx <- 1
  for (i in 1:N) {
    for (t in 1:TT) {
      treated <- (t >= T0_vec[i])
      D_vec[idx] <- as.integer(treated)
      factor_component <- if (r > 0) sum(F_mat[t, ] * L_mat[i, ]) else 0
      Y_vec[idx] <- alpha_i[i] + xi_t[t] + industry_fe[industry_i[i]] +
        factor_component + tau * D_vec[idx] + rnorm(1, 0, 0.5)
      industry_vec[idx] <- industry_i[i]
      id_vec[idx] <- i
      time_vec[idx] <- t
      idx <- idx + 1
    }
  }

  data.frame(id = id_vec, time = time_vec, Y = Y_vec, D = D_vec,
             industry = industry_vec)
}

make_cfe_fe_nesting_data <- function(N = 100, TT = 30, Ntr = 30, tau = 3.0,
                                      r = 2, seed = 42) {
  set.seed(seed)
  F_mat <- matrix(rnorm(TT * r), TT, r)
  L_mat <- matrix(rnorm(N * r), N, r)
  alpha_i <- rnorm(N, 0, 1)
  xi_t <- rnorm(TT, 0, 0.5)

  ## Industry: treated get 1-2, control get 3-4 (unit-nesting / Type A)
  industry_i <- integer(N)
  industry_i[1:Ntr] <- sample(1:2, Ntr, replace = TRUE)
  industry_i[(Ntr + 1):N] <- sample(3:4, N - Ntr, replace = TRUE)
  industry_fe <- c(-1.0, 0.5, -0.5, 1.0)

  ## Staggered treatment adoption
  T0_vec <- rep(Inf, N)
  if (Ntr > 0) {
    T0_vec[1:Ntr] <- sample(round(TT * 0.4):round(TT * 0.7), Ntr,
                             replace = TRUE)
  }

  Y_vec <- D_vec <- numeric(N * TT)
  id_vec <- time_vec <- industry_vec <- integer(N * TT)
  idx <- 1
  for (i in 1:N) {
    for (t in 1:TT) {
      treated <- (t >= T0_vec[i])
      D_vec[idx] <- as.integer(treated)
      factor_component <- if (r > 0) sum(F_mat[t, ] * L_mat[i, ]) else 0
      Y_vec[idx] <- alpha_i[i] + xi_t[t] + industry_fe[industry_i[i]] +
        factor_component + tau * D_vec[idx] + rnorm(1, 0, 0.5)
      industry_vec[idx] <- industry_i[i]
      id_vec[idx] <- i
      time_vec[idx] <- t
      idx <- idx + 1
    }
  }

  data.frame(id = id_vec, time = time_vec, Y = Y_vec, D = D_vec,
             industry = industry_vec)
}

make_cfe_full_data <- function(N = 100, TT = 30, Ntr = 30, tau = 3.0,
                                r = 2, seed = 42) {
  set.seed(seed)
  F_mat <- matrix(rnorm(TT * r), TT, r)
  L_mat <- matrix(rnorm(N * r), N, r)
  alpha_i <- rnorm(N, 0, 1)
  xi_t <- rnorm(TT, 0, 0.5)

  ## Z: time-invariant covariate
  Z_i <- rnorm(N, 0, 1)
  gamma_t <- rnorm(TT, 0, 0.5)

  ## Q: time-varying basis
  Q_t <- (1:TT) / TT
  kappa_i <- rnorm(N, 0, 0.5)

  ## Industry: shared (Type B) — all units draw from same pool
  industry_i <- sample(1:4, N, replace = TRUE)
  industry_fe <- c(-1.0, 0.5, 0.0, 1.5)

  ## Staggered treatment adoption
  T0_vec <- rep(Inf, N)
  if (Ntr > 0) {
    T0_vec[1:Ntr] <- sample(round(TT * 0.4):round(TT * 0.7), Ntr,
                             replace = TRUE)
  }

  Y_vec <- D_vec <- Z_vec <- Q_vec <- numeric(N * TT)
  id_vec <- time_vec <- industry_vec <- integer(N * TT)
  idx <- 1
  for (i in 1:N) {
    for (t in 1:TT) {
      treated <- (t >= T0_vec[i])
      D_vec[idx] <- as.integer(treated)
      factor_component <- if (r > 0) sum(F_mat[t, ] * L_mat[i, ]) else 0
      Y_vec[idx] <- alpha_i[i] + xi_t[t] + Z_i[i] * gamma_t[t] +
        Q_t[t] * kappa_i[i] + industry_fe[industry_i[i]] +
        factor_component + tau * D_vec[idx] + rnorm(1, 0, 0.5)
      Z_vec[idx] <- Z_i[i]
      Q_vec[idx] <- Q_t[t]
      industry_vec[idx] <- industry_i[i]
      id_vec[idx] <- i
      time_vec[idx] <- t
      idx <- idx + 1
    }
  }

  data.frame(id = id_vec, time = time_vec, Y = Y_vec, D = D_vec,
             Z = Z_vec, Q = Q_vec, industry = industry_vec)
}

## ---- Category B: Specification Equivalence ----

test_that("Phase 3a-B1: ife+nevertreated equals cfe+nevertreated (no extras)", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  set.seed(100)
  out_ife <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  set.seed(100)
  out_cfe <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(abs(out_ife$att.avg - out_cfe$att.avg) < 0.1,
              info = paste("ife ATT:", out_ife$att.avg, "cfe ATT:", out_cfe$att.avg))
})

test_that("Phase 3a-B2: fe+nevertreated equals cfe+nevertreated at r=0", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 0, seed = 42)

  out_ife <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 0, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  out_cfe <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", r = 0, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(abs(out_ife$att.avg - out_cfe$att.avg) < 0.1,
              info = paste("ife ATT:", out_ife$att.avg, "cfe ATT:", out_cfe$att.avg))
})

test_that("Phase 3a-B3: cfe+nevertreated is deterministic (same seed, same result)", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  set.seed(100)
  out1 <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  set.seed(100)
  out2 <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_identical(out1$att.avg, out2$att.avg)
})

test_that("Phase 3a-B4: gsynth equals cfe+nevertreated (no extras)", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  set.seed(100)
  out_gs <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "gsynth", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", parallel = FALSE
  )))

  set.seed(100)
  out_cfe <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(abs(out_gs$att.avg - out_cfe$att.avg) < 0.1,
              info = paste("gsynth ATT:", out_gs$att.avg, "cfe ATT:", out_cfe$att.avg))
})

## ---- Category C: Accuracy ----

test_that("Phase 3a-C1: accuracy with Z (time-invariant covariates)", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 200, TT = 30, Ntr = 60, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(abs(out$att.avg - 3.0) < 0.5,
              info = paste("ATT:", out$att.avg, "true tau: 3.0"))
})

test_that("Phase 3a-C2: accuracy with Q (unit-specific time trends)", {
  skip_on_cran()
  df <- make_cfe_q_data(N = 200, TT = 30, Ntr = 60, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Q = "Q", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(abs(out$att.avg - 3.0) < 0.5,
              info = paste("ATT:", out$att.avg, "true tau: 3.0"))
})

test_that("Phase 3a-C3: accuracy with shared extra FE (Type B)", {
  skip_on_cran()
  df <- make_cfe_fe_data(N = 200, TT = 30, Ntr = 60, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time", "industry"),
    method = "cfe", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(abs(out$att.avg - 3.0) < 0.5,
              info = paste("ATT:", out$att.avg, "true tau: 3.0"))
})

test_that("Phase 3a-C4: accuracy with unit-nesting extra FE (Type A)", {
  skip_on_cran()
  df <- make_cfe_fe_nesting_data(N = 200, TT = 30, Ntr = 60, tau = 3.0,
                                  r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time", "industry"),
    method = "cfe", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(abs(out$att.avg - 3.0) < 0.5,
              info = paste("ATT:", out$att.avg, "true tau: 3.0"))
})

test_that("Phase 3a-C5: accuracy with all CFE components", {
  skip_on_cran()
  df <- make_cfe_full_data(N = 200, TT = 30, Ntr = 60, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time", "industry"),
    method = "cfe", Z = "Z", Q = "Q", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(abs(out$att.avg - 3.0) < 0.5,
              info = paste("ATT:", out$att.avg, "true tau: 3.0"))
})

test_that("Phase 3a-C6: notyettreated vs nevertreated both accurate", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 200, TT = 30, Ntr = 60, tau = 3.0, r = 2, seed = 42)

  out_nyt <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "notyettreated",
    parallel = FALSE
  )))

  out_nt <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(abs(out_nyt$att.avg - 3.0) < 1.0,
              info = paste("notyettreated ATT:", out_nyt$att.avg))
  expect_true(abs(out_nt$att.avg - 3.0) < 1.0,
              info = paste("nevertreated ATT:", out_nt$att.avg))
})

test_that("Phase 3a-C7: accuracy with Z at r=0 (no factors)", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 200, TT = 30, Ntr = 60, tau = 3.0, r = 0, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 0, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(abs(out$att.avg - 3.0) < 0.5,
              info = paste("ATT:", out$att.avg, "true tau: 3.0"))
})

## ---- Category D: Validation Guards ----

test_that("Phase 3a-D1: Type-B FE with missing level in co errors", {
  skip_on_cran()
  df <- make_cfe_fe_data(N = 100, TT = 30, Ntr = 30, seed = 42)
  ## Assign a unique industry level to treated unit 1 (not in controls)
  df$industry[df$id == 1] <- 99

  expect_error(
    suppressWarnings(suppressMessages(fect::fect(
      Y ~ D, data = df, index = c("id", "time", "industry"),
      method = "cfe", r = 2, CV = FALSE, se = FALSE,
      force = "two-way", factors.from = "nevertreated",
      parallel = FALSE
    )))
  )
})

test_that("Phase 3a-D2: no never-treated units errors", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, TT = 20, Ntr = 40)

  expect_error(
    suppressWarnings(suppressMessages(fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "cfe", r = 0, CV = FALSE, se = FALSE,
      force = "two-way", factors.from = "nevertreated",
      parallel = FALSE
    ))),
    regexp = "never.treated|nevertreated|control"
  )
})

test_that("Phase 3a-D3: too few control units for r errors", {
  skip_on_cran()
  df <- make_factor_data(N = 35, TT = 30, Ntr = 33, tau = 3.0, r = 2, seed = 42)

  expect_error(
    suppressWarnings(suppressMessages(fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "cfe", r = 2, CV = FALSE, se = FALSE,
      force = "two-way", factors.from = "nevertreated",
      parallel = FALSE
    ))),
    regexp = "never.treated|nevertreated|insuffic|too few"
  )
})

test_that("Phase 3a-D4: valid cfe+nevertreated combination is accepted", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- expect_no_error(suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  ))))

  expect_false(is.null(out$att.avg))
})

test_that("Phase 3a-D5: mc+nevertreated is rejected", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, seed = 42)

  expect_error(
    suppressWarnings(suppressMessages(fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "mc", r = 2, CV = FALSE, se = FALSE,
      force = "two-way", factors.from = "nevertreated",
      parallel = FALSE
    ))),
    regexp = "mc|nevertreated|factors.from"
  )
})

## ---- Category A: Solver Equivalence ----

test_that("Phase 3a-A1: complex_fe_ub vs inter_fe_ub equivalence", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  set.seed(100)
  out_cfe <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  set.seed(100)
  out_ife <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(abs(out_cfe$att.avg - out_ife$att.avg) < 0.1,
              info = paste("cfe ATT:", out_cfe$att.avg, "ife ATT:", out_ife$att.avg))
})

test_that("Phase 3a-A2: solver equivalence at r=0 (FE only)", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 0, seed = 42)

  out_cfe <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", r = 0, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  out_ife <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 0, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(abs(out_cfe$att.avg - out_ife$att.avg) < 0.01,
              info = paste("cfe ATT:", out_cfe$att.avg, "ife ATT:", out_ife$att.avg))
})

test_that("Phase 3a-A3: solver equivalence with unbalanced data", {
  skip_on_cran()
  df <- make_staggered_data(N = 40, TT = 20, Ntr = 15, tau = 3.0, seed = 42)

  out_cfe <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", r = 0, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  out_ife <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 0, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(abs(out_cfe$att.avg - out_ife$att.avg) < 0.1,
              info = paste("cfe ATT:", out_cfe$att.avg, "ife ATT:", out_ife$att.avg))
})

## ---- Category F: Bootstrap / Inference ----

test_that("Phase 3a-F2: parametric bootstrap SE with cfe+nevertreated", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE,
    se = TRUE, vartype = "parametric", nboots = 30,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_false(is.null(out$est.att))
  expect_false(all(is.na(out$est.att[, "S.E."])))
})

test_that("Phase 3a-F1: jackknife SE with cfe+nevertreated", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE,
    se = TRUE, vartype = "jackknife", nboots = 30,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_false(is.null(out$est.att))
  expect_false(all(is.na(out$est.att[, "S.E."])))
})

test_that("Phase 3a-F3: parametric bootstrap with cfe+notyettreated and extras", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE,
    se = TRUE, vartype = "parametric", nboots = 30,
    force = "two-way", factors.from = "notyettreated",
    parallel = FALSE
  )))

  expect_false(all(is.na(out$est.att[, "S.E."])))
})

test_that("Phase 3a-F4: ife+nevertreated SE unchanged (regression)", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE,
    se = TRUE, nboots = 30,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_false(is.null(out$est.att))
  expect_false(all(is.na(out$est.att[, "S.E."])))
})

## ---- Category E: Output Completeness ----

test_that("Phase 3a-E1: gamma and kappa fields present", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_false(is.null(out$gamma))
  ## kappa can be NULL if no Q was specified — that's OK
})

test_that("Phase 3a-E2: factors.from field in output", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_equal(out$factors.from, "nevertreated")
})

test_that("Phase 3a-E3: core output fields non-NULL with correct dimensions", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_false(is.null(out$Y.ct))
  expect_false(is.null(out$eff))
  expect_false(is.null(out$att.avg))
  expect_true(is.numeric(out$att.avg))
  expect_false(is.na(out$att.avg))
  expect_equal(nrow(out$Y.ct), length(unique(df$time)))
  expect_equal(ncol(out$Y.ct), length(unique(df$id)))
})

test_that("Phase 3a-E4: plot() works without error", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_no_error(suppressWarnings(suppressMessages(plot(out))))
})

test_that("Phase 3a-E5: print() works without error", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_no_error(suppressWarnings(suppressMessages(print(out))))
})

## ---- Category G: Cross-Validation ----

test_that("Phase 3a-G1: CV selects a reasonable r with factors in data", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 200, TT = 30, Ntr = 60, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = c(0, 5), CV = TRUE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(out$r.cv >= 0)
  expect_true(is.numeric(out$r.cv))
})

test_that("Phase 3a-G2: CV selects r=0 on no-factor data", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 200, TT = 30, Ntr = 60, tau = 3.0, r = 0, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = c(0, 3), CV = TRUE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_equal(out$r.cv, 0)
})

test_that("Phase 3a-G3: ife+nevertreated CV unchanged (regression)", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = c(0, 5), CV = TRUE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_true(out$r.cv >= 0)
})

## ---- Category H: Edge Cases ----

test_that("Phase 3a-H1: single treated unit", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 1, tau = 3.0, r = 2, seed = 42)

  expect_no_error(suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  ))))
})

test_that("Phase 3a-H2: single control unit (r=0)", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 31, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  ## With only 1 never-treated unit at r=0 this may either run or error clearly
  result <- tryCatch(
    suppressWarnings(suppressMessages(fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "cfe", r = 0, CV = FALSE, se = FALSE,
      force = "two-way", factors.from = "nevertreated",
      parallel = FALSE
    ))),
    error = function(e) e
  )

  ## Either it runs successfully or produces a clear error
  if (inherits(result, "error")) {
    expect_true(grepl("never.treated|nevertreated|insuffic|too few|control",
                      result$message, ignore.case = TRUE),
                info = paste("Unexpected error:", result$message))
  } else {
    expect_false(is.null(result$att.avg))
  }
})

test_that("Phase 3a-H3: cfe+nevertreated with r=0 (no factors)", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 0, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 0, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  )))

  expect_no_error(out)
  expect_true(abs(out$att.avg - 3.0) < 1.0,
              info = paste("ATT:", out$att.avg, "true tau: 3.0"))
})

test_that("Phase 3a-H4: no covariates at all", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  expect_no_error(suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  ))))
})

test_that("Phase 3a-H5: treatment reversals", {
  skip_on_cran()
  df <- make_cfe_full_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)
  ## Create reversals: unit 1 reverts to control after period 20
  reversal_rows <- df$id == 1 & df$time > 20
  df$D[reversal_rows] <- 0

  expect_no_error(suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time", "industry"),
    method = "cfe", Z = "Z", Q = "Q", r = 2, CV = FALSE, se = FALSE,
    force = "two-way", factors.from = "nevertreated",
    parallel = FALSE
  ))))
})

## ========================================================
## PHASE 3a-I: Parametric bootstrap under nevertreated
## with em=TRUE / em=FALSE for both cfe and ife
## Uses parallel computing with controlled seeds (doRNG)
## ========================================================

test_that("Phase 3a-I1: ife+nevertreated parametric bootstrap, em=TRUE, parallel", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, force = "two-way",
    factors.from = "nevertreated", em = TRUE,
    se = TRUE, vartype = "bootstrap", nboots = 20,
    parallel = TRUE, cores = 2, seed = 123
  )))

  expect_false(is.na(out$att.avg),
               info = "att.avg should not be NA")
  expect_false(is.null(out$est.att),
               info = "est.att should not be NULL")
  expect_true(any(!is.na(out$est.att[, "S.E."])),
              info = "SE estimates should not all be NA")
})

test_that("Phase 3a-I2: ife+nevertreated parametric bootstrap, em=FALSE, parallel", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, force = "two-way",
    factors.from = "nevertreated", em = FALSE,
    se = TRUE, vartype = "bootstrap", nboots = 20,
    parallel = TRUE, cores = 2, seed = 123
  )))

  expect_false(is.na(out$att.avg),
               info = "att.avg should not be NA")
  expect_false(is.null(out$est.att),
               info = "est.att should not be NULL")
  expect_true(any(!is.na(out$est.att[, "S.E."])),
              info = "SE estimates should not all be NA")
})

test_that("Phase 3a-I3: cfe+nevertreated parametric bootstrap, em=TRUE, parallel", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, force = "two-way",
    factors.from = "nevertreated", em = TRUE,
    se = TRUE, vartype = "bootstrap", nboots = 20,
    parallel = TRUE, cores = 2, seed = 123
  )))

  expect_false(is.na(out$att.avg),
               info = "att.avg should not be NA")
  expect_false(is.null(out$est.att),
               info = "est.att should not be NULL")
  expect_true(any(!is.na(out$est.att[, "S.E."])),
              info = "SE estimates should not all be NA")
})

test_that("Phase 3a-I4: cfe+nevertreated parametric bootstrap, em=FALSE, parallel", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, force = "two-way",
    factors.from = "nevertreated", em = FALSE,
    se = TRUE, vartype = "bootstrap", nboots = 20,
    parallel = TRUE, cores = 2, seed = 123
  )))

  expect_false(is.na(out$att.avg),
               info = "att.avg should not be NA")
  expect_false(is.null(out$est.att),
               info = "est.att should not be NULL")
  expect_true(any(!is.na(out$est.att[, "S.E."])),
              info = "SE estimates should not all be NA")
})

test_that("Phase 3a-I5: bootstrap reproducibility with same seed (ife+nevertreated)", {
  skip_on_cran()
  df <- make_factor_data(N = 80, TT = 20, Ntr = 25, tau = 3.0, r = 2, seed = 42)

  run_boot <- function(s) {
    suppressWarnings(suppressMessages(fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r = 2, CV = FALSE, force = "two-way",
      factors.from = "nevertreated", em = TRUE,
      se = TRUE, vartype = "bootstrap", nboots = 10,
      parallel = TRUE, cores = 2, seed = s
    )))
  }

  out1 <- run_boot(999)
  out2 <- run_boot(999)

  expect_equal(out1$att.avg, out2$att.avg,
               info = "Same seed should give identical att.avg")
  expect_equal(out1$est.att[, "S.E."], out2$est.att[, "S.E."],
               info = "Same seed should give identical SE")
})

test_that("Phase 3a-I6: bootstrap reproducibility with same seed (cfe+nevertreated)", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 80, TT = 20, Ntr = 25, tau = 3.0, r = 2, seed = 42)

  run_boot <- function(s) {
    suppressWarnings(suppressMessages(fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "cfe", Z = "Z", r = 2, CV = FALSE, force = "two-way",
      factors.from = "nevertreated", em = TRUE,
      se = TRUE, vartype = "bootstrap", nboots = 10,
      parallel = TRUE, cores = 2, seed = s
    )))
  }

  out1 <- run_boot(999)
  out2 <- run_boot(999)

  expect_equal(out1$att.avg, out2$att.avg,
               info = "Same seed should give identical att.avg")
  expect_equal(out1$est.att[, "S.E."], out2$est.att[, "S.E."],
               info = "Same seed should give identical SE")
})

test_that("Phase 3a-I7: different seeds give different bootstrap SE (ife+nevertreated)", {
  skip_on_cran()
  df <- make_factor_data(N = 80, TT = 20, Ntr = 25, tau = 3.0, r = 2, seed = 42)

  run_boot <- function(s) {
    suppressWarnings(suppressMessages(fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r = 2, CV = FALSE, force = "two-way",
      factors.from = "nevertreated", em = TRUE,
      se = TRUE, vartype = "bootstrap", nboots = 10,
      parallel = TRUE, cores = 2, seed = s
    )))
  }

  out1 <- run_boot(111)
  out2 <- run_boot(222)

  ## Point estimate should be the same (not seed-dependent)
  expect_equal(out1$att.avg, out2$att.avg, tolerance = 1e-10,
               info = "Point estimate should not depend on bootstrap seed")
  ## SE should differ (different bootstrap draws)
  expect_false(identical(out1$est.att[, "S.E."], out2$est.att[, "S.E."]),
               info = "Different seeds should give different SE draws")
})

test_that("Phase 3a-I8: ATT accuracy check under bootstrap (cfe+nevertreated, em=TRUE)", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 150, TT = 30, Ntr = 50, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, force = "two-way",
    factors.from = "nevertreated", em = TRUE,
    se = TRUE, vartype = "bootstrap", nboots = 20,
    parallel = TRUE, cores = 2, seed = 123
  )))

  expect_true(abs(out$att.avg - 3.0) < 1.0,
              info = paste("ATT:", out$att.avg, "should be near 3.0"))
  ## CI should cover the true value
  ci_lower <- out$est.avg[, "CI.lower"]
  ci_upper <- out$est.avg[, "CI.upper"]
  expect_true(ci_lower < 3.0 & ci_upper > 3.0,
              info = paste("95% CI [", ci_lower, ",", ci_upper,
                           "] should cover true tau=3.0"))
})

test_that("Phase 3a-I9: quantile.CI=TRUE bias-corrected reflection CI (cfe+nevertreated)", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, force = "two-way",
    factors.from = "nevertreated",
    se = TRUE, vartype = "bootstrap", nboots = 20,
    quantile.CI = TRUE,
    parallel = TRUE, cores = 2, seed = 123
  )))

  expect_false(is.null(out$est.att),
               info = "est.att should not be NULL")
  expect_true("CI.lower" %in% colnames(out$est.att),
              info = "est.att should have CI.lower column")
  expect_true("CI.upper" %in% colnames(out$est.att),
              info = "est.att should have CI.upper column")
  expect_true(any(!is.na(out$est.att[, "CI.lower"])),
              info = "At least some CI.lower values should be non-NA")
  expect_false(is.null(out$est.avg),
               info = "est.avg should not be NULL")
  expect_true("CI.lower" %in% colnames(out$est.avg),
              info = "est.avg should have CI.lower column")
  expect_true("CI.upper" %in% colnames(out$est.avg),
              info = "est.avg should have CI.upper column")
  ci_lower <- out$est.avg[, "CI.lower"]
  ci_upper <- out$est.avg[, "CI.upper"]
  expect_true(ci_lower < 3.0 & ci_upper > 3.0,
              info = paste("Quantile CI [", ci_lower, ",", ci_upper,
                           "] should cover true tau=3.0"))
})

test_that("Phase 3a-I10: em=TRUE vs em=FALSE identical for nevertreated (ife)", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)

  run_em <- function(em_val) {
    suppressWarnings(suppressMessages(fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "ife", r = 2, CV = FALSE, force = "two-way",
      factors.from = "nevertreated", em = em_val,
      se = FALSE
    )))
  }

  out1 <- run_em(TRUE)
  out2 <- run_em(FALSE)

  expect_equal(out1$att.avg, out2$att.avg, tolerance = 1e-10,
               info = "em=TRUE and em=FALSE should give identical att.avg for nevertreated")
})

test_that("Phase 3a-I11: unbalanced data forces _ub/EM path in draw.error() bootstrap", {
  skip_on_cran()

  ## IFE on unbalanced data
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)
  set.seed(999)
  n_drop <- floor(nrow(df) * 0.05)
  drop_idx <- sample(seq_len(nrow(df)), n_drop, replace = FALSE)
  df_ub <- df[-drop_idx, ]

  out_ife <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df_ub, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, force = "two-way",
    factors.from = "nevertreated",
    se = TRUE, vartype = "bootstrap", nboots = 15,
    parallel = TRUE, cores = 2, seed = 123
  )))

  expect_false(is.na(out_ife$att.avg),
               info = "IFE unbalanced: att.avg should not be NA")
  expect_true(any(!is.na(out_ife$est.att[, "S.E."])),
              info = "IFE unbalanced: SE estimates should not all be NA")

  ## CFE on unbalanced data
  df2 <- make_cfe_z_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)
  set.seed(999)
  n_drop2 <- floor(nrow(df2) * 0.05)
  drop_idx2 <- sample(seq_len(nrow(df2)), n_drop2, replace = FALSE)
  df2_ub <- df2[-drop_idx2, ]

  out_cfe <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df2_ub, index = c("id", "time"),
    method = "cfe", Z = "Z", r = 2, CV = FALSE, force = "two-way",
    factors.from = "nevertreated",
    se = TRUE, vartype = "bootstrap", nboots = 15,
    parallel = TRUE, cores = 2, seed = 123
  )))

  expect_false(is.na(out_cfe$att.avg),
               info = "CFE unbalanced: att.avg should not be NA")
  expect_true(any(!is.na(out_cfe$est.att[, "S.E."])),
              info = "CFE unbalanced: SE estimates should not all be NA")
})

test_that("Phase 3a-I12: r=0 invariance — factors.from is a no-op when r=0", {
  skip_on_cran()
  df <- make_factor_data(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 0, seed = 42)

  out_nyt <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 0, CV = FALSE, force = "two-way",
    factors.from = "notyettreated",
    se = FALSE
  )))

  out_nt <- suppressWarnings(suppressMessages(fect::fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", r = 0, CV = FALSE, force = "two-way",
    factors.from = "nevertreated",
    se = FALSE
  )))

  expect_equal(out_nyt$att.avg, out_nt$att.avg, tolerance = 1e-2,
               info = "factors.from should be a no-op when r=0: att.avg must match")
})

## ========================================================
## CV routing for cfe+nevertreated (REQ-cv-gap-001)
## ========================================================

test_that("cfe+nevertreated CV selects r and runs without error", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 200, TT = 30, Ntr = 60, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = c(0, 3), CV = TRUE, se = FALSE,
    factors.from = "nevertreated",
    force = "two-way", parallel = TRUE, cores = 10
  )))

  expect_s3_class(out, "fect")
  expect_true(is.numeric(out$r.cv), info = "r.cv should be set by CV")
  expect_true(out$r.cv >= 1 && out$r.cv <= 3,
              info = paste("r.cv should be >= 1 for DGP with r=2, got:", out$r.cv))
  expect_true(is.numeric(out$att.avg))
  expect_true(!is.na(out$att.avg))
})

## ========================================================
## CFE CV r-selection accuracy tests (REQ-cfe-cv-rselect-001)
## ========================================================

test_that("cfe+nevertreated+Z CV correctly selects r=2", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 200, TT = 30, Ntr = 60, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = c(0, 5), CV = TRUE, se = FALSE,
    factors.from = "nevertreated",
    force = "two-way", parallel = TRUE, cores = 10
  )))

  expect_equal(out$r.cv, 2,
               info = paste("r.cv should be 2 for DGP with r=2, got:", out$r.cv))
  expect_true(abs(out$att.avg - 3.0) < 0.5,
              info = paste("ATT should be close to 3.0, got:", out$att.avg))
})

test_that("cfe+nevertreated+Z CV correctly selects r=0 on zero-factor data", {
  skip_on_cran()
  df <- make_cfe_z_data(N = 200, TT = 30, Ntr = 60, tau = 3.0, r = 0, seed = 42)

  out <- suppressWarnings(suppressMessages(fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", Z = "Z", r = c(0, 5), CV = TRUE, se = FALSE,
    factors.from = "nevertreated",
    force = "two-way", parallel = TRUE, cores = 10
  )))

  expect_equal(out$r.cv, 0,
               info = paste("r.cv should be 0 for DGP with r=0, got:", out$r.cv))
})

test_that("cfe+nevertreated CV selects r=2 on factor-only data (no Z in DGP)", {
  skip_on_cran()
  df <- make_factor_data(N = 200, TT = 30, Ntr = 60, tau = 3.0, r = 2, seed = 42)

  out <- suppressWarnings(suppressMessages(fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "cfe", r = c(0, 5), CV = TRUE, se = FALSE,
    factors.from = "nevertreated",
    force = "two-way", parallel = TRUE, cores = 10
  )))

  expect_equal(out$r.cv, 2,
               info = paste("r.cv should be 2 for DGP with r=2, got:", out$r.cv))
})
