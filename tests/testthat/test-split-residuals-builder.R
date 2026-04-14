## Unit tests for split_residuals parameter — builder verification
## Tests based on spec.md §3–§14 (Stage 4 debiasing POC).
## Do NOT reference test-spec.md — builder tests are independent.

## NOTE: these tests are for builder self-verification only.
## Tester will run its own independent test suite (test-paraboot-splitting.R).

suppressMessages(library(fect))

## --------------------------------------------------------------------------
## Shared helper: small balanced panel with Nco >= 10
## --------------------------------------------------------------------------
make_panel <- function(N = 30, TT = 20, Ntr = 6, seed = 101) {
  set.seed(seed)
  id   <- rep(1:N, each = TT)
  time <- rep(1:TT, N)
  D    <- rep(0L, N * TT)
  for (u in seq_len(Ntr)) {
    idx  <- which(id == u & time >= 11)
    D[idx] <- 1L
  }
  Y <- rnorm(N * TT) + D * 2
  data.frame(id = id, time = time, Y = Y, D = D)
}

df30 <- make_panel(N = 30, TT = 20, Ntr = 6)

## --------------------------------------------------------------------------
## Test 1: partition_controls() returns correct structure
## --------------------------------------------------------------------------
test_that("partition_controls returns two disjoint halves covering id.co", {
  set.seed(7)
  id.co <- c(10L, 15L, 20L, 25L, 30L, 35L, 40L, 45L)  # Nco = 8
  res <- fect:::partition_controls(id.co, K = 2L)

  expect_type(res, "list")
  expect_named(res, c("A", "B"))

  # Union covers all of id.co
  expect_setequal(c(res$A, res$B), id.co)

  # Disjoint
  expect_length(intersect(res$A, res$B), 0)

  # Sizes: floor(8/2) = 4 for A, 4 for B
  expect_length(res$A, 4L)
  expect_length(res$B, 4L)

  # Each half is sorted
  expect_identical(res$A, sort(res$A))
  expect_identical(res$B, sort(res$B))
})

## --------------------------------------------------------------------------
## Test 2: partition_controls() — odd Nco: A is smaller, B is larger
## --------------------------------------------------------------------------
test_that("partition_controls with odd Nco gives floor/ceiling split", {
  set.seed(13)
  id.co <- 1:9  # Nco = 9 (odd)
  res <- fect:::partition_controls(id.co, K = 2L)

  expect_length(res$A, 4L)  # floor(9/2) = 4
  expect_length(res$B, 5L)  # ceil(9/2)  = 5
  expect_setequal(c(res$A, res$B), id.co)
})

## --------------------------------------------------------------------------
## Test 3: partition_controls() — Nco < 4 raises error
## --------------------------------------------------------------------------
test_that("partition_controls stops when Nco < 4", {
  expect_error(
    fect:::partition_controls(1:3, K = 2L),
    regexp = "split_residuals requires at least 4 control units"
  )
  expect_error(
    fect:::partition_controls(1:1, K = 2L),
    regexp = "split_residuals requires at least 4 control units"
  )
})

## --------------------------------------------------------------------------
## Test 4: partition_controls() — K != 2 raises error
## --------------------------------------------------------------------------
test_that("partition_controls stops for K != 2", {
  expect_error(
    fect:::partition_controls(1:8, K = 3L),
    regexp = "only K=2 is supported"
  )
})

## --------------------------------------------------------------------------
## Test 5: split_residuals validation — bad type raises error
## --------------------------------------------------------------------------
test_that("split_residuals must be logical scalar", {
  skip_if_not_installed("fect")
  expect_error(
    fect(Y ~ D, data = df30, index = c("id", "time"), method = "ife",
         se = TRUE, vartype = "parametric",
         time.component.from = "nevertreated",
         nboots = 10, parallel = FALSE, seed = 1,
         split_residuals = "yes"),
    regexp = "split_residuals.*must be TRUE or FALSE"
  )
  expect_error(
    fect(Y ~ D, data = df30, index = c("id", "time"), method = "ife",
         se = TRUE, vartype = "parametric",
         time.component.from = "nevertreated",
         nboots = 10, parallel = FALSE, seed = 1,
         split_residuals = NA),
    regexp = "split_residuals.*must be TRUE or FALSE"
  )
})

## --------------------------------------------------------------------------
## Test 6: split_residuals=FALSE (explicit) — parity with no argument
## --------------------------------------------------------------------------
test_that("split_residuals=FALSE produces identical output to not specifying it", {
  r_default <- fect(Y ~ D, data = df30, index = c("id", "time"),
                    method = "ife", se = TRUE, vartype = "parametric",
                    time.component.from = "nevertreated",
                    nboots = 50, parallel = FALSE, seed = 999)

  r_false   <- fect(Y ~ D, data = df30, index = c("id", "time"),
                    method = "ife", se = TRUE, vartype = "parametric",
                    time.component.from = "nevertreated",
                    nboots = 50, parallel = FALSE, seed = 999,
                    split_residuals = FALSE)

  expect_identical(r_default$att.avg,  r_false$att.avg)
  expect_identical(r_default$att.avg.W, r_false$att.avg.W)
  expect_identical(r_default$att,       r_false$att)
})

## --------------------------------------------------------------------------
## Test 7: Gate C fires with split_residuals=FALSE + notyettreated
## --------------------------------------------------------------------------
test_that("Gate C fires when split_residuals=FALSE and notyettreated", {
  expect_error(
    fect(Y ~ D, data = df30, index = c("id", "time"), method = "ife",
         se = TRUE, vartype = "parametric",
         time.component.from = "notyettreated",
         nboots = 10, parallel = FALSE, seed = 1,
         split_residuals = FALSE),
    regexp = "debiased parametric bootstrap"
  )
})

## --------------------------------------------------------------------------
## Test 8: Gate C bypassed with split_residuals=TRUE + notyettreated;
##         experimental message is emitted
## --------------------------------------------------------------------------
test_that("split_residuals=TRUE bypasses Gate C and emits experimental message", {
  msgs <- character(0)
  res <- withCallingHandlers(
    fect(Y ~ D, data = df30, index = c("id", "time"), method = "ife",
         se = TRUE, vartype = "parametric",
         time.component.from = "notyettreated",
         nboots = 30, parallel = FALSE, seed = 42,
         split_residuals = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  # No error
  expect_false(inherits(res, "error"))
  # Experimental message was emitted
  expect_true(any(grepl("experimental", msgs, ignore.case = TRUE)))
  # Valid output structure
  expect_true(!is.null(res$att.avg))
})

## --------------------------------------------------------------------------
## Test 9: split_residuals=TRUE + nevertreated runs without error
## --------------------------------------------------------------------------
test_that("split_residuals=TRUE works with nevertreated path", {
  expect_no_error(
    fect(Y ~ D, data = df30, index = c("id", "time"), method = "ife",
         se = TRUE, vartype = "parametric",
         time.component.from = "nevertreated",
         nboots = 30, parallel = FALSE, seed = 42,
         split_residuals = TRUE)
  )
})

## --------------------------------------------------------------------------
## Test 10: split_residuals=TRUE with Nco < 4 raises partition error
## --------------------------------------------------------------------------
test_that("split_residuals=TRUE with Nco=3 raises partition_controls error", {
  # Panel with only 3 controls
  df_tiny <- make_panel(N = 6, TT = 15, Ntr = 3)
  expect_error(
    fect(Y ~ D, data = df_tiny, index = c("id", "time"), method = "ife",
         se = TRUE, vartype = "parametric",
         time.component.from = "notyettreated",
         nboots = 5, parallel = FALSE, seed = 1,
         split_residuals = TRUE),
    regexp = "split_residuals requires at least 4 control units"
  )
})

## --------------------------------------------------------------------------
## Test 11: split_residuals=TRUE produces wider CIs than split=FALSE
##          (direction test — not a hard assertion, just a smoke direction check)
## --------------------------------------------------------------------------
test_that("split_residuals=TRUE produces valid CI output (structure check)", {
  r_split <- fect(Y ~ D, data = df30, index = c("id", "time"),
                  method = "ife", se = TRUE, vartype = "parametric",
                  time.component.from = "nevertreated",
                  nboots = 50, parallel = FALSE, seed = 42,
                  split_residuals = TRUE)

  # CI.lower and CI.upper exist and are finite
  expect_true(!is.null(r_split$est.avg))
  ci_lower <- r_split$est.avg[, "CI.lower"]
  ci_upper <- r_split$est.avg[, "CI.upper"]
  expect_true(all(is.finite(ci_lower) | is.na(ci_lower)))
  expect_true(all(is.finite(ci_upper) | is.na(ci_upper)))
  # CIs open upward: lower < upper wherever both are finite
  finite_idx <- is.finite(ci_lower) & is.finite(ci_upper)
  if (any(finite_idx)) {
    expect_true(all(ci_lower[finite_idx] < ci_upper[finite_idx]))
  }
})
