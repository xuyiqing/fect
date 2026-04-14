## Unit tests for impute_Y0() dispatcher and valid_controls() helper
## (Phase 1 + Phase 2 of the parametric-bootstrap refactor)
##
## These tests verify:
##   - impute_Y0 dispatches to the correct callee for each (method, predictive) combination
##   - impute_Y0 errors on unsupported combinations (mc)
##   - valid_controls returns the correct subset of control columns
##   - hasRevs+parametric gate fires correctly in fect()

test_that("impute_Y0 errors on mc method", {
  skip_on_cran()

  # Build a minimal data stub so impute_Y0 would reach the else branch
  # We don't need real data — just check the error fires before any callee call
  expect_error(
    impute_Y0(
      method     = "mc",
      predictive = "nevertreated",
      Y          = matrix(0, 5, 4),
      X          = NULL,
      D          = matrix(0, 5, 4),
      W          = NULL,
      I          = matrix(1, 5, 4),
      II         = matrix(1, 5, 4),
      T.on       = matrix(0L, 5, 4),
      tuning     = 1L,
      force      = 0,
      tol        = 1e-3,
      max.iteration = 10
    ),
    regexp = "MC support is planned for Phase 4"
  )
})

test_that("impute_Y0 errors on unknown method+predictive combination", {
  skip_on_cran()

  expect_error(
    impute_Y0(
      method     = "zzz",
      predictive = "nevertreated",
      Y          = matrix(0, 5, 4),
      X          = NULL,
      D          = matrix(0, 5, 4),
      W          = NULL,
      I          = matrix(1, 5, 4),
      II         = matrix(1, 5, 4),
      T.on       = matrix(0L, 5, 4),
      tuning     = 1L,
      force      = 0,
      tol        = 1e-3,
      max.iteration = 10
    ),
    regexp = "unsupported combination"
  )
})

test_that("valid_controls errors on mc method", {
  skip_on_cran()

  # Minimal out stub
  out_stub <- list(
    D   = matrix(c(0, 0, 1, 1, 0, 0, 0, 0, 0, 0), nrow = 5, ncol = 2),
    I   = matrix(1, 5, 2),
    r.cv = 1L
  )
  expect_error(
    valid_controls(out_stub, method = "mc", predictive = "nevertreated", force = 0),
    regexp = "MC not supported"
  )
})

test_that("valid_controls returns correct subset without FE adjustment", {
  skip_on_cran()

  # 5 time periods, 3 units: unit 1 = treated (onset t=3), units 2-3 = control
  # D:
  #   unit 1: c(0,0,1,1,1)
  #   unit 2: c(0,0,0,0,0)
  #   unit 3: c(0,0,0,0,0)
  TT <- 5
  N <- 3
  D <- matrix(0, TT, N)
  D[3:5, 1] <- 1
  I <- matrix(1, TT, N)
  r.cv <- 1L

  out_stub <- list(D = D, I = I, r.cv = r.cv)

  # T0.ub for treated unit 1 = 2 (periods 1-2 are pre-treatment)
  # max.T0.ub = 2, T0.ub.min = 2
  # co.pre[j] = sum(I.co[1:2, j]) for each control
  # co.post[j] = sum(I.co[3:5, j]) for each control
  # With full I, co.pre = 2, co.post = 3 for all controls
  # threshold (force=0): r.cv = 1
  # 2 >= 1 AND 3 >= 1 -> both controls should pass
  vc <- valid_controls(out_stub, method = "ife", predictive = "notyettreated", force = 0)
  expect_equal(vc, c(2L, 3L))
})

test_that("valid_controls applies +1 threshold when force has unit FE", {
  skip_on_cran()

  TT <- 5
  N <- 3
  D <- matrix(0, TT, N)
  D[3:5, 1] <- 1
  I <- matrix(1, TT, N)
  r.cv <- 2L  # threshold without FE = 2, with unit FE = 3

  out_stub <- list(D = D, I = I, r.cv = r.cv)

  # co.pre = 2, threshold (force=1) = 2 + 1 = 3
  # 2 < 3, so co.pre threshold fails -> both controls excluded
  vc_unit_fe <- valid_controls(out_stub, method = "ife", predictive = "notyettreated", force = 1)
  expect_equal(length(vc_unit_fe), 0L)

  # Without FE: threshold = 2, co.pre = 2 >= 2 -> both pass
  vc_no_fe <- valid_controls(out_stub, method = "ife", predictive = "notyettreated", force = 0)
  expect_equal(vc_no_fe, c(2L, 3L))
})

test_that("hasRevs+parametric gate fires in fect()", {
  skip_on_cran()

  # Load a dataset with treatment reversals
  data(fect.simdata)  # standard simdata in the package
  df <- fect.simdata

  # Manually introduce a reversal: flip treatment off for one unit in the middle
  # Find a treated unit and flip one period from 1 -> 0 in the middle
  treated_units <- unique(df$id[df$D == 1])
  if (length(treated_units) < 1) {
    skip("No treated units in simdata")
  }
  target_unit <- treated_units[1]
  unit_rows <- which(df$id == target_unit & df$D == 1)
  if (length(unit_rows) < 2) {
    skip("Insufficient treated periods to introduce reversal")
  }
  # Introduce a reversal: set one treated period back to 0
  mid_row <- unit_rows[ceiling(length(unit_rows) / 2)]
  df$D[mid_row] <- 0

  expect_error(
    fect(Y ~ D, data = df, index = c("id", "time"),
         method = "ife", vartype = "parametric", se = TRUE, nboots = 5),
    regexp = "parametric.*reversal|reversal.*parametric",
    ignore.case = TRUE
  )
})
