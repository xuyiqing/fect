## ---------------------------------------------------------------
## Tests for the $sample logical mask on the fect() return value.
##
## $sample is a logical matrix (same dims as $Y.dat) that is TRUE for
## every cell used in any part of the estimation procedure:
##   obs.missing == 1  (treated)
##   obs.missing == 2  (control / pre-treatment)
##   obs.missing == 5  (placebo or carryover period)
## Cells with obs.missing == 3 (missing/unbalanced) or 4 (removed unit)
## are FALSE.
##
## Derived via:  sample <- obs.missing %in% c(1L, 2L, 5L)
## ---------------------------------------------------------------

suppressWarnings(data("simdata", package = "fect"))

## -- S.1  Slot exists and is a logical matrix --

test_that("S.1: fit$sample exists and is a logical matrix", {

  skip_on_cran()

  set.seed(1)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
      method = "fe", CV = FALSE, se = FALSE, parallel = FALSE
    )
  ))

  expect_true("sample" %in% names(fit))
  expect_true(is.matrix(fit$sample))
  expect_true(is.logical(fit$sample))
})


## -- S.2  Dimensions match Y.dat --

test_that("S.2: fit$sample has the same dimensions as fit$Y.dat", {

  skip_on_cran()

  set.seed(1)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
      method = "fe", CV = FALSE, se = FALSE, parallel = FALSE
    )
  ))

  expect_equal(dim(fit$sample), dim(fit$Y.dat))
})


## -- S.3  Consistency with obs.missing codes --
## sample must equal obs.missing %in% c(1L, 2L, 5L) cell-for-cell.

test_that("S.3: fit$sample equals obs.missing %in% c(1L, 2L, 5L)", {

  skip_on_cran()

  set.seed(1)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
      method = "fe", CV = FALSE, se = FALSE, parallel = FALSE
    )
  ))

  expected <- matrix(fit$obs.missing %in% c(1L, 2L, 5L),
                     nrow = nrow(fit$obs.missing), ncol = ncol(fit$obs.missing),
                     dimnames = dimnames(fit$obs.missing))
  expect_identical(fit$sample, expected)
})


## -- S.4  Removed units are FALSE --
## Units with rm.id (obs.missing == 4) must be entirely FALSE in $sample.

test_that("S.4: cells with obs.missing == 4 (removed) are FALSE in $sample", {

  skip_on_cran()

  set.seed(1)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
      method = "fe", CV = FALSE, se = FALSE, parallel = FALSE
    )
  ))

  removed_cells <- fit$obs.missing == 4L
  ## If no removed units, skip the body (not an error — just no removed units).
  if (any(removed_cells)) {
    expect_true(all(!fit$sample[removed_cells]))
  } else {
    expect_true(TRUE)  ## no removed units — pass trivially
  }
})


## -- S.5  Missing cells are FALSE --
## obs.missing == 3 cells must be FALSE in $sample.

test_that("S.5: cells with obs.missing == 3 (missing) are FALSE in $sample", {

  skip_on_cran()

  set.seed(1)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
      method = "fe", CV = FALSE, se = FALSE, parallel = FALSE
    )
  ))

  missing_cells <- fit$obs.missing == 3L
  if (any(missing_cells)) {
    expect_true(all(!fit$sample[missing_cells]))
  } else {
    expect_true(TRUE)
  }
})


## -- S.6  Placebo period cells are TRUE --
## When placeboTest = TRUE, pre-treatment placebo cells (obs.missing == 5)
## must appear as TRUE in $sample.

test_that("S.6: placebo cells (obs.missing == 5) are TRUE in $sample", {

  skip_on_cran()

  set.seed(2)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
      method = "fe", CV = FALSE, se = FALSE, parallel = FALSE,
      placeboTest = TRUE, placebo.period = c(-2, -1)
    )
  ))

  placebo_cells <- fit$obs.missing == 5L
  if (any(placebo_cells)) {
    expect_true(all(fit$sample[placebo_cells]))
  } else {
    ## If no placebo cells found, the placeboTest may not have produced code 5
    ## under this data/seed — at minimum check slot still exists and is logical.
    expect_true(is.logical(fit$sample))
  }
})


## -- S.7  at_least_some_TRUE: basic sanity that $sample is not all FALSE --

test_that("S.7: fit$sample has at least some TRUE cells (treated + control)", {

  skip_on_cran()

  set.seed(1)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
      method = "fe", CV = FALSE, se = FALSE, parallel = FALSE
    )
  ))

  expect_true(any(fit$sample))
})


## -- S.8  Treated cells are TRUE --
## obs.missing == 1 cells must be TRUE in $sample.

test_that("S.8: treated cells (obs.missing == 1) are TRUE in $sample", {

  skip_on_cran()

  set.seed(1)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
      method = "fe", CV = FALSE, se = FALSE, parallel = FALSE
    )
  ))

  treated_cells <- fit$obs.missing == 1L
  expect_true(any(treated_cells))
  expect_true(all(fit$sample[treated_cells]))
})


## -- S.9  Control cells are TRUE --
## obs.missing == 2 cells must be TRUE in $sample.

test_that("S.9: control cells (obs.missing == 2) are TRUE in $sample", {

  skip_on_cran()

  set.seed(1)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
      method = "fe", CV = FALSE, se = FALSE, parallel = FALSE
    )
  ))

  control_cells <- fit$obs.missing == 2L
  expect_true(any(control_cells))
  expect_true(all(fit$sample[control_cells]))
})


## -- S.10  Method = ife gives same slot structure --

test_that("S.10: fit$sample present and logical for method='ife'", {

  skip_on_cran()

  set.seed(3)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
      method = "ife", r = 2, CV = FALSE, se = FALSE, parallel = FALSE
    )
  ))

  expect_true("sample" %in% names(fit))
  expect_true(is.logical(fit$sample))
  expect_equal(dim(fit$sample), dim(fit$Y.dat))
})
