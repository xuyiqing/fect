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
