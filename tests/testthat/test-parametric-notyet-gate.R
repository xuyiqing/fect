## test-parametric-notyet-gate.R
## Gate tests for parametric + notyettreated hard gate
## StatsClaw run: 2026-04-14-notyet-parametric-gate
##
## Positive-fire tests (NGP-1 through NGP-4): gate fires when it should.
## Negative-fire tests (NGN-1 through NGN-3): gate does NOT fire for valid paths.

## ---- helpers ----------------------------------------------------------------

make_gate_fixture <- function(seed = 42) {
  ## Staggered adoption panel with 50 units (10 treated, 40 never-treated)
  ## Uses same structure as make_fixture_data() in test-paraboot-parity.R
  set.seed(seed)
  N <- 50; TT <- 20; T0 <- 12; Ntr <- 10
  id   <- rep(1:N, each = TT)
  time <- rep(1:TT, N)
  D    <- as.integer(id <= Ntr & time >= T0)
  alpha_i <- rep(rnorm(N, 0, 1), each = TT)
  xi_t    <- rep(rnorm(TT, 0, 0.5), N)
  Y <- alpha_i + xi_t + 2 * D + rnorm(N * TT, 0, 0.5)
  data.frame(id = id, time = time, Y = Y, D = D)
}

## Verify fixture has never-treated units (required for nevertreated tests)
.check_nev <- function(d) {
  D_mat <- tapply(d$D, list(d$time, d$id), identity)
  id_co <- which(colSums(D_mat) == 0)
  stopifnot("Fixture must have never-treated units" = length(id_co) > 0)
  invisible(TRUE)
}

## ============================================================
## POSITIVE-FIRE TESTS: gate fires for invalid combinations
## ============================================================

test_that("NGP-1: ife + vartype='parametric' + default factors.from errors", {
  skip_on_cran()
  d <- make_gate_fixture()
  expect_error(
    fect(
      Y ~ D, data = d, index = c("id", "time"),
      method = "ife", r = 1, se = TRUE,
      vartype = "parametric", nboots = 10,
      CV = FALSE, parallel = FALSE
    ),
    regexp = "parametric",
    ignore.case = TRUE
  )
})

test_that("NGP-2: cfe + vartype='parametric' + default factors.from errors", {
  skip_on_cran()
  d <- make_gate_fixture()
  expect_error(
    fect(
      Y ~ D, data = d, index = c("id", "time"),
      method = "cfe", r = 1, se = TRUE,
      vartype = "parametric", nboots = 10,
      CV = FALSE, parallel = FALSE
    ),
    regexp = "parametric",
    ignore.case = TRUE
  )
})

test_that("NGP-3: ife + vartype='parametric' + explicit factors.from='notyettreated' errors", {
  skip_on_cran()
  d <- make_gate_fixture()
  expect_error(
    fect(
      Y ~ D, data = d, index = c("id", "time"),
      method = "ife", time.component.from = "notyettreated",
      r = 1, se = TRUE, vartype = "parametric", nboots = 10,
      CV = FALSE, parallel = FALSE
    ),
    regexp = "parametric",
    ignore.case = TRUE
  )
})

test_that("NGP-4: error message mentions both an alternative (nevertreated or bootstrap or jackknife)", {
  skip_on_cran()
  d <- make_gate_fixture()
  expect_error(
    fect(
      Y ~ D, data = d, index = c("id", "time"),
      method = "ife", r = 1, se = TRUE,
      vartype = "parametric", nboots = 10,
      CV = FALSE, parallel = FALSE
    ),
    regexp = "nevertreated|bootstrap|jackknife",
    ignore.case = TRUE
  )
})

## ============================================================
## NEGATIVE-FIRE TESTS: gate does NOT fire for valid paths
## ============================================================

test_that("NGN-1: gsynth + parametric does NOT error", {
  skip_on_cran()
  d <- make_gate_fixture()
  out <- suppressWarnings(suppressMessages(fect(
    Y ~ D, data = d, index = c("id", "time"),
    method = "gsynth", r = 1, se = TRUE,
    vartype = "parametric", nboots = 10,
    CV = FALSE, parallel = FALSE
  )))
  expect_s3_class(out, "fect")
})

test_that("NGN-2: ife + nevertreated + parametric does NOT error", {
  skip_on_cran()
  d <- make_gate_fixture()
  .check_nev(d)
  out <- suppressWarnings(suppressMessages(fect(
    Y ~ D, data = d, index = c("id", "time"),
    method = "ife", time.component.from = "nevertreated",
    r = 1, se = TRUE, vartype = "parametric", nboots = 10,
    CV = FALSE, parallel = FALSE
  )))
  expect_s3_class(out, "fect")
})

test_that("NGN-3: cfe + nevertreated + parametric does NOT error", {
  skip_on_cran()
  d <- make_gate_fixture()
  .check_nev(d)
  out <- suppressWarnings(suppressMessages(fect(
    Y ~ D, data = d, index = c("id", "time"),
    method = "cfe", time.component.from = "nevertreated",
    r = 1, se = TRUE, vartype = "parametric", nboots = 10,
    CV = FALSE, parallel = FALSE
  )))
  expect_s3_class(out, "fect")
})
