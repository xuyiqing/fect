## ---------------------------------------------------------------
## Tests for imputed_outcomes() and the slot contract.
##
## Coverage: column schema, cells filter (NULL / logical / formula),
## replicates expansion, direction = "on" / "off", error paths.
## ---------------------------------------------------------------

suppressWarnings(data("simdata", package = "fect"))


.fit_canonical <- function(keep.sims = FALSE, nboots = 50) {
  set.seed(42)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = simdata, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, nboots = nboots, parallel = FALSE,
      keep.sims = keep.sims
    )
  ))
}


## -- I.1  Column schema, point estimates ----------------------------

test_that("I.1: imputed_outcomes() returns the documented column schema", {

  skip_on_cran()

  fit <- .fit_canonical()
  po  <- fect::imputed_outcomes(fit)

  expect_s3_class(po, "data.frame")

  expected_cols <- c("id", "time", "event.time", "cohort", "treated",
                     "Y_obs", "Y0_hat", "eff", "eff_debias", "W.agg")
  expect_setequal(names(po), expected_cols)

  ## Every row is treated
  expect_true(all(po$treated))

  ## Y0_hat = Y_obs - eff (within numeric tolerance)
  expect_equal(po$Y0_hat, po$Y_obs - po$eff)

  ## eff_debias is zero for plain imputation estimators (FE)
  expect_true(all(po$eff_debias == 0))

  ## W.agg is 1 when fit has no W
  expect_true(all(po$W.agg == 1))

  ## event.time is populated at most treated cells; can be NA for
  ## units that begin in the treated state (no pre-treatment data to
  ## anchor the event-time grid). Document and accept.
  expect_true(mean(!is.na(po$event.time)) > 0.95)
})


## -- I.2  cells filter: NULL = unfiltered ---------------------------

test_that("I.2: cells = NULL returns all treated cells", {

  skip_on_cran()

  fit <- .fit_canonical()
  po_full  <- fect::imputed_outcomes(fit, cells = NULL)
  po_again <- fect::imputed_outcomes(fit)
  expect_identical(po_full, po_again)
})


## -- I.3  cells filter: logical vector ------------------------------

test_that("I.3: cells = logical vector subsets correctly", {

  skip_on_cran()

  fit <- .fit_canonical()
  po  <- fect::imputed_outcomes(fit)

  mask <- po$event.time %in% 1:5
  po_logical <- fect::imputed_outcomes(fit, cells = mask)

  expect_equal(nrow(po_logical), sum(mask))
  expect_true(all(po_logical$event.time %in% 1:5))
})


## -- I.4  cells filter: one-sided formula ---------------------------

test_that("I.4: cells = formula subsets correctly", {

  skip_on_cran()

  fit <- .fit_canonical()
  po_form  <- fect::imputed_outcomes(fit, cells = ~ event.time %in% 1:5)
  po_full  <- fect::imputed_outcomes(fit)

  expect_equal(nrow(po_form), sum(po_full$event.time %in% 1:5))
  expect_true(all(po_form$event.time %in% 1:5))
})


## -- I.5  replicates = TRUE returns each replicate's own cells ------
## (2.4.6: replicate b's rows are the copies of the treated cells that the
## resampled panel of replicate b contains, read through colnames.boot. Under
## the case bootstrap a unit drawn twice appears twice and an undrawn unit
## not at all, so the row count varies by replicate.)

test_that("I.5: replicates = TRUE returns each replicate's cells with a replicate column", {

  skip_on_cran()

  fit <- .fit_canonical(keep.sims = TRUE, nboots = 25)
  po  <- fect::imputed_outcomes(fit)
  po_b <- fect::imputed_outcomes(fit, replicates = TRUE)

  expect_true("replicate" %in% names(po_b))
  expect_setequal(unique(po_b$replicate), seq_len(25))
  expect_identical(names(po_b), c(names(po), "replicate"))

  ## Averaging eff within a replicate gives that replicate's att.avg
  ## (unweighted fit, all treated cells).
  rep_mean <- c(tapply(po_b$eff, po_b$replicate, mean))
  expect_equal(unname(rep_mean), unname(c(fit$att.avg.boot)),
               tolerance = 1e-10)

  ## Every row is a treated cell of the point surface, and Y0_hat = Y_obs - eff
  expect_true(all(paste(po_b$id, po_b$time) %in% paste(po$id, po$time)))
  expect_equal(po_b$Y0_hat, po_b$Y_obs - po_b$eff)

  ## Copies of one cell share (id, time, Y_obs, W.agg); Y0_hat varies
  ## across replicates.
  key <- paste(po_b$id, po_b$time)
  busiest <- names(which.max(table(key)))
  one_cell <- po_b[key == busiest, ]
  expect_gt(length(unique(one_cell$replicate)), 1)
  expect_equal(length(unique(one_cell$Y_obs)), 1)
  expect_equal(length(unique(one_cell$W.agg)), 1)
  expect_gt(length(unique(one_cell$Y0_hat)), 1)
})


## -- I.6  replicates = TRUE without keep.sims errors with locked wording

test_that("I.6: replicates = TRUE without keep.sims errors with locked wording", {

  skip_on_cran()

  fit <- .fit_canonical(keep.sims = FALSE)
  expect_error(
    fect::imputed_outcomes(fit, replicates = TRUE),
    "keep\\.sims = TRUE",
    fixed = FALSE
  )
})


## -- I.7  direction = "off" works on reversal panels ----------------

test_that("I.7: direction = 'off' returns a sensible long form on reversal panels", {

  skip_on_cran()

  fit <- .fit_canonical()
  ## simdata has reversal; fit$T.off should be populated.
  expect_false(is.null(fit$T.off))

  po_off <- fect::imputed_outcomes(fit, direction = "off")
  expect_true(nrow(po_off) > 0)
  expect_true(all(!is.na(po_off$event.time)))

  ## Cohort under "off" is the calendar time of first exit.
  expect_true(all(!is.na(po_off$cohort) | po_off$event.time != 1))
})


## -- I.8  validate_po_contract: bad fit errors ----------------------

test_that("I.8: passing a non-fect object errors with a contract message", {

  expect_error(
    fect::imputed_outcomes("not a fit"),
    "fect"
  )
})
