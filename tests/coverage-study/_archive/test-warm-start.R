## Tests for v2.4.3 partial warm-start API.
##
## Design: warm.start = c("none", "linear"). "linear" enables partial warm-start
## in the bootstrap loop (auxiliaries warmed from main fit, factors cold-start).
##
## Acceptance criteria from ref/v242-warm-start-investigation/partial-warm-design.md:
## - point estimates byte-identical between cold and warm
## - SE diff < 5% relative when both at tol = 1e-5
## - tol > 1e-5 + warm.start = "linear" rejected with informative error

skip_if_not_installed("fect")
data(simdata, package = "fect")

## ---- W1. Default is "none" (no behavior change for existing scripts) ----
test_that("W1: warm.start default is 'none' (no regression for pre-2.4.3 scripts)", {
    set.seed(1)
    fit_default <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                        method = "ife", r = 2, force = "two-way",
                        se = FALSE, CV = FALSE)
    set.seed(1)
    fit_explicit <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                         method = "ife", r = 2, force = "two-way",
                         se = FALSE, CV = FALSE,
                         warm.start = "none")
    expect_identical(fit_default$att.avg, fit_explicit$att.avg)
    expect_identical(fit_default$est.att[, "ATT"], fit_explicit$est.att[, "ATT"])
})

## ---- W2. Invalid value rejected ----
test_that("W2: invalid warm.start value is rejected", {
    expect_error(
        fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
             method = "ife", r = 2, force = "two-way",
             se = FALSE, CV = FALSE,
             warm.start = "full"),
        regexp = "must be one of"
    )
    expect_error(
        fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
             method = "ife", r = 2, force = "two-way",
             se = FALSE, CV = FALSE,
             warm.start = c("none","linear")),
        regexp = "must be one of"
    )
})

## ---- W3. tol > 1e-5 + warm.start = "linear" rejected ----
## Note: v2.4.3 default tol is now 1e-5, so we must explicitly pass a
## looser tol to trigger the validation.
test_that("W3: warm.start = 'linear' requires tol <= 1e-5", {
    expect_error(
        fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
             method = "ife", r = 2, force = "two-way",
             se = TRUE, vartype = "bootstrap", nboots = 5,
             parallel = FALSE, CV = FALSE,
             tol = 1e-3,                  ## explicit loose tol
             warm.start = "linear"),
        regexp = "warm\\.start = 'linear' requires tol <= 1e-5"
    )
})

## ---- W4. warm.start = "linear" runs without error at tol = 1e-5 (IFE) ----
test_that("W4: warm.start = 'linear' runs end-to-end on IFE at tol = 1e-5", {
    skip_on_cran()
    set.seed(1)
    fit <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                method = "ife", r = 2, force = "two-way",
                se = TRUE, vartype = "bootstrap", nboots = 10,
                parallel = FALSE, CV = FALSE,
                tol = 1e-5, max.iteration = 5000,
                warm.start = "linear")
    expect_true(is.numeric(fit$att.avg))
    expect_true(is.finite(fit$att.avg))
    expect_true(all(c("ATT", "S.E.") %in% colnames(fit$est.att)))
})

## ---- W5. Point estimate byte-identical to cold (warm.start affects bootstrap only) ----
test_that("W5: warm.start = 'linear' preserves point estimate", {
    skip_on_cran()
    set.seed(1)
    fit_cold <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                     method = "ife", r = 2, force = "two-way",
                     se = TRUE, vartype = "bootstrap", nboots = 10,
                     parallel = FALSE, CV = FALSE,
                     tol = 1e-5, max.iteration = 5000)
    set.seed(1)
    fit_warm <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                     method = "ife", r = 2, force = "two-way",
                     se = TRUE, vartype = "bootstrap", nboots = 10,
                     parallel = FALSE, CV = FALSE,
                     tol = 1e-5, max.iteration = 5000,
                     warm.start = "linear")
    expect_equal(fit_cold$att.avg, fit_warm$att.avg, tolerance = 1e-10)
    expect_equal(fit_cold$est.att[, "ATT"],
                 fit_warm$est.att[, "ATT"], tolerance = 1e-10)
})

## ---- W6. CFE smoke test ----
test_that("W6: warm.start = 'linear' runs on CFE method", {
    skip_on_cran()
    set.seed(1)
    fit <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
                method = "cfe", r = 2, force = "two-way",
                se = TRUE, vartype = "bootstrap", nboots = 10,
                parallel = FALSE, CV = FALSE,
                tol = 1e-5, max.iteration = 5000,
                warm.start = "linear")
    expect_true(is.numeric(fit$att.avg))
    expect_true(is.finite(fit$att.avg))
})
