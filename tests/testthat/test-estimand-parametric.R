## -------------------------------------------------------------------
## v2.4.1 adds "parametric" to the vartype enum on estimand(). The
## boot.R machinery already populates fit$eff.boot under
## vartype = "parametric" + keep.sims = TRUE; the v2.4.1 change is just
## the match.arg gate. These tests verify:
##   - the new enum value is accepted
##   - byte-equality with fit$est.att under parametric
##   - all four `type` paths (att, att.cumu, aptt, log.att) work on
##     parametric fits, sourcing from fit$eff.boot
##   - keep.sims = FALSE still errors helpfully on non-fast paths
##   - vartype = "none" returns NA SE/CI under parametric
## -------------------------------------------------------------------

## Helper: cache one parametric fit per session for the heavy tests.
.make_parametric_fit <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    skip_on_cran()
    e <- new.env()
    data(sim_linear, package = "fect", envir = e)
    set.seed(42)
    fit <- fect(Y ~ D, data = e$sim_linear, index = c("id", "time"),
                method = "ife", force = "two-way", se = TRUE,
                nboots = 30, r = 2, CV = FALSE, keep.sims = TRUE,
                vartype = "parametric",
                time.component.from = "nevertreated",
                parallel = FALSE)
    cached <<- fit
    fit
  }
})

test_that("vartype = 'BAD' is rejected (enum gate, CRAN-runnable)", {
  ## Cheap arg-validation. match.arg fails at the top of estimand()
  ## before any fit-shape inspection, so a stub list is sufficient.
  stub <- structure(list(), class = "fect")
  expect_error(
    estimand(stub, "att", "event.time", vartype = "BAD"),
    "should be one of"
  )
})

test_that("vartype = 'parametric' is accepted by estimand()", {
  skip_on_cran()
  fit <- .make_parametric_fit()
  res <- estimand(fit, "att", "event.time", vartype = "parametric")
  expect_s3_class(res, "data.frame")
  expect_true(all(res$vartype == "parametric"))
})

test_that("estimand(fit, 'att', 'event.time') byte-equals fit$est.att under parametric", {
  skip_on_cran()
  fit <- .make_parametric_fit()
  res <- estimand(fit, "att", "event.time")
  expect_identical(res$estimate, unname(fit$est.att[, "ATT"]))
  expect_identical(res$se,       unname(fit$est.att[, "S.E."]))
  expect_identical(res$ci.lo,    unname(fit$est.att[, "CI.lower"]))
  expect_identical(res$ci.hi,    unname(fit$est.att[, "CI.upper"]))
  expect_identical(res$n_cells,  unname(fit$est.att[, "count"]))
})

test_that("output 'vartype' column reports fit$vartype, not the argument", {
  skip_on_cran()
  fit <- .make_parametric_fit()
  ## Default arg defaults to "bootstrap"; output should still say
  ## "parametric" because that's what generated the replicates.
  res_default <- estimand(fit, "att", "event.time")
  expect_true(all(res_default$vartype == "parametric"))
  ## Explicit "parametric": same numerical result, same vartype column.
  res_explicit <- estimand(fit, "att", "event.time", vartype = "parametric")
  expect_identical(res_default, res_explicit)
})

test_that("att.cumu / aptt / att-overall work on parametric fits", {
  skip_on_cran()
  fit <- .make_parametric_fit()

  expect_silent(r1 <- estimand(fit, "att.cumu", "event.time"))
  expect_true(all(r1$vartype == "parametric"))
  expect_true(any(!is.na(r1$estimate)))

  expect_silent(r2 <- estimand(fit, "att", "overall", window = c(1, 5)))
  expect_true(r2$vartype == "parametric")
  expect_true(!is.na(r2$estimate))

  expect_silent(r3 <- estimand(fit, "aptt", "event.time"))
  expect_true(all(r3$vartype == "parametric"))
  expect_true(any(!is.na(r3$estimate)))
})

test_that("log.att hard-errors on parametric fits with negative Y (v2.4.2+)", {
  skip_on_cran()
  fit <- .make_parametric_fit()
  ## sim_linear has many negative Y cells, so log.att now hard-errors
  ## on the bootstrap cell-drop pathology (v2.4.2+). The previous
  ## v2.4.1 behavior of silently warning + dropping cells produced
  ## meaningless inference; the hard-error redirects users to the
  ## actionable options.
  expect_error(
    estimand(fit, "log.att", "event.time"),
    "log-ATT bootstrap is unreliable"
  )
})

test_that("vartype = 'none' under parametric returns NA SE/CI", {
  skip_on_cran()
  fit <- .make_parametric_fit()
  res <- estimand(fit, "att", "overall", window = c(1, 5),
                  vartype = "none")
  expect_true(!is.na(res$estimate))
  expect_true(is.na(res$se))
  expect_true(is.na(res$ci.lo))
  expect_true(is.na(res$ci.hi))
})

test_that("keep.sims = FALSE under parametric still errors helpfully on non-fast paths", {
  skip_on_cran()
  e <- new.env()
  data(sim_linear, package = "fect", envir = e)
  set.seed(42)
  fit_no_sims <- fect(Y ~ D, data = e$sim_linear, index = c("id", "time"),
                      method = "ife", force = "two-way", se = TRUE,
                      nboots = 20, r = 2, CV = FALSE,
                      keep.sims = FALSE,
                      vartype = "parametric",
                      time.component.from = "nevertreated",
                      parallel = FALSE)
  expect_null(fit_no_sims$eff.boot)

  ## att / event.time fast path reads fit$est.att, so still works.
  expect_silent(res_fast <- estimand(fit_no_sims, "att", "event.time"))
  expect_true(all(res_fast$vartype == "parametric"))

  ## att.cumu requires eff.boot; should error with locked wording.
  expect_error(
    estimand(fit_no_sims, "att.cumu", "event.time"),
    "No bootstrap/jackknife results available"
  )
})

test_that("eff.boot replicate count integrity under parametric", {
  skip_on_cran()
  fit <- .make_parametric_fit()
  expect_equal(dim(fit$eff.boot)[3], 30L)
  expect_equal(dim(fit$eff.boot)[1], nrow(fit$Y.dat))
  expect_equal(dim(fit$eff.boot)[2], ncol(fit$Y.dat))
})
