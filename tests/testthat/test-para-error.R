## Unit tests for para.error argument (spec: wild-as-paraerror)
##
## Tests verify:
##   1. para.error = "auto" default produces correct resolved mode stored on fit
##   2. para.error = "empirical" runs without error on fully-observed panel
##   3. para.error = "wild" runs without error on fully-observed panel, all 5
##      ci.methods return CIs containing the true ATT
##   4. vartype = "wild" (deprecated alias) emits a warning and produces
##      fit$vartype == "parametric" and fit$para.error == "wild"
##   5. para.error = "wild" hard-errors on missing-data panel
##   6. para.error = "empirical" hard-errors on missing-data panel
##   7. Invalid para.error value produces a clear error message
##   8. para.error is ignored when vartype != "parametric"
##   9. para.error = "ar" works on fully-observed panel

## -----------------------------------------------------------------------
## Shared test fixture
## -----------------------------------------------------------------------

make_panel <- function(N = 40, TT = 20, T0 = 12, Ntr = 12, seed = 101,
                       add_missing = FALSE, n_missing = 50) {
    set.seed(seed)
    alpha_i <- rnorm(N, 0, 2)
    xi_t    <- rnorm(TT, 0, 1)
    D       <- matrix(0L, TT, N)
    D[(T0 + 1):TT, 1:Ntr] <- 1L
    eps <- matrix(rnorm(N * TT, 0, 1), TT, N)
    Y   <- outer(xi_t, rep(1, N)) + outer(rep(1, TT), alpha_i) + 3.0 * D + eps
    d   <- data.frame(
        id   = rep(1:N, each = TT),
        time = rep(1:TT, N),
        Y    = as.vector(Y),
        D    = as.vector(D)
    )
    if (add_missing) {
        set.seed(seed + 999)
        d$Y[sample(nrow(d), n_missing)] <- NA
    }
    d
}

fect_para <- function(d, para.error = "auto", vartype = "parametric",
                      nboots = 100, seed = 42, ...) {
    set.seed(seed)
    suppressWarnings(suppressMessages(
        fect(
            Y ~ D,
            data   = d,
            index  = c("id", "time"),
            method = "fe",
            force  = "two-way",
            time.component.from = "nevertreated",
            se      = TRUE,
            vartype = vartype,
            para.error = para.error,
            nboots     = nboots,
            parallel   = FALSE,
            keep.sims  = TRUE,
            CV         = FALSE,
            ...
        )
    ))
}

## -----------------------------------------------------------------------
## Test 1: auto default resolves to "empirical" for fully-observed panel
## -----------------------------------------------------------------------

test_that("para.error = 'auto' resolves to 'empirical' on fully-observed panel", {
    skip_on_cran()
    d   <- make_panel()
    fit <- fect_para(d, para.error = "auto")
    expect_equal(fit$para.error, "empirical")
    expect_equal(fit$vartype, "parametric")
})

## -----------------------------------------------------------------------
## Test 2: para.error = "empirical" runs and stores correct label
## -----------------------------------------------------------------------

test_that("para.error = 'empirical' runs on fully-observed panel", {
    skip_on_cran()
    d   <- make_panel()
    fit <- fect_para(d, para.error = "empirical")
    expect_equal(fit$para.error, "empirical")

    ## CI should be sensible (contains true ATT ≈ 3)
    res <- estimand(fit, "att", "overall", window = c(1, 8), ci.method = "bca")
    expect_true(res$ci.lo < 3.0 && res$ci.hi > 3.0,
        info = sprintf("BCA CI [%.3f, %.3f] should contain 3.0", res$ci.lo, res$ci.hi))
})

## -----------------------------------------------------------------------
## Test 3: para.error = "wild" runs and all ci.methods return sensible CIs
## -----------------------------------------------------------------------

test_that("para.error = 'wild' runs on fully-observed panel and all ci.methods work", {
    skip_on_cran()
    d   <- make_panel()
    fit <- fect_para(d, para.error = "wild", nboots = 200)
    expect_equal(fit$para.error, "wild")
    expect_equal(fit$vartype, "parametric")

    for (m in c("normal", "basic", "percentile", "bc", "bca")) {
        res <- estimand(fit, "att", "overall", window = c(1, 8), ci.method = m)
        expect_true(
            !is.na(res$ci.lo) && !is.na(res$ci.hi),
            info = sprintf("ci.method='%s': CI should not be NA", m)
        )
        expect_true(
            res$ci.lo < 3.0 && res$ci.hi > 3.0,
            info = sprintf(
                "ci.method='%s': CI [%.3f, %.3f] should contain 3.0",
                m, res$ci.lo, res$ci.hi
            )
        )
        expect_true(
            res$ci.hi - res$ci.lo > 0,
            info = sprintf("ci.method='%s': CI width should be positive", m)
        )
    }
})

## -----------------------------------------------------------------------
## Test 4: vartype = "wild" deprecated alias — warning emitted, correct routing
## -----------------------------------------------------------------------

test_that("vartype = 'wild' emits deprecation warning and routes to parametric + wild", {
    skip_on_cran()
    d <- make_panel()
    fit <- NULL
    expect_warning(
        {
            set.seed(42)
            fit <- suppressMessages(fect(
                Y ~ D,
                data   = d,
                index  = c("id", "time"),
                method = "fe",
                force  = "two-way",
                time.component.from = "nevertreated",
                se      = TRUE,
                vartype = "wild",
                nboots  = 100,
                parallel = FALSE,
                keep.sims = TRUE,
                CV = FALSE
            ))
        },
        regexp = "deprecated",
        fixed  = FALSE
    )
    expect_equal(fit$vartype,    "parametric")
    expect_equal(fit$para.error, "wild")
})

## -----------------------------------------------------------------------
## Test 5: para.error = "wild" hard-errors on missing-data panel
## -----------------------------------------------------------------------

test_that("para.error = 'wild' hard-errors on missing-data panel", {
    skip_on_cran()
    d <- make_panel(add_missing = TRUE)
    expect_error(
        fect(
            Y ~ D,
            data   = d,
            index  = c("id", "time"),
            method = "fe",
            force  = "two-way",
            time.component.from = "nevertreated",
            se      = TRUE,
            vartype = "parametric",
            para.error = "wild",
            na.rm   = FALSE,
            nboots  = 50,
            parallel = FALSE,
            CV = FALSE
        ),
        regexp = "fully-observed panel",
        fixed  = FALSE
    )
})

## -----------------------------------------------------------------------
## Test 6: para.error = "empirical" hard-errors on missing-data panel
## -----------------------------------------------------------------------

test_that("para.error = 'empirical' hard-errors on missing-data panel", {
    skip_on_cran()
    d <- make_panel(add_missing = TRUE)
    expect_error(
        fect(
            Y ~ D,
            data   = d,
            index  = c("id", "time"),
            method = "fe",
            force  = "two-way",
            time.component.from = "nevertreated",
            se      = TRUE,
            vartype = "parametric",
            para.error = "empirical",
            na.rm   = FALSE,
            nboots  = 50,
            parallel = FALSE,
            CV = FALSE
        ),
        regexp = "fully-observed panel",
        fixed  = FALSE
    )
})

## -----------------------------------------------------------------------
## Test 7: Invalid para.error value produces clear error
## -----------------------------------------------------------------------

test_that("invalid para.error value produces clear error", {
    skip_on_cran()
    d <- make_panel()
    expect_error(
        fect(
            Y ~ D,
            data   = d,
            index  = c("id", "time"),
            method = "fe",
            force  = "two-way",
            time.component.from = "nevertreated",
            se      = TRUE,
            vartype = "parametric",
            para.error = "invalid_value",
            nboots  = 50,
            parallel = FALSE,
            CV = FALSE
        ),
        regexp = "para.error",
        fixed  = FALSE
    )
})

## -----------------------------------------------------------------------
## Test 8: para.error is ignored when vartype != "parametric"
## -----------------------------------------------------------------------

test_that("para.error is accepted but ignored for vartype = 'bootstrap'", {
    skip_on_cran()
    d   <- make_panel()
    ## Should run without error even though para.error = "wild" is passed
    ## with vartype = "bootstrap" (non-parametric path ignores it)
    fit <- NULL
    expect_no_error({
        set.seed(42)
        fit <- suppressWarnings(suppressMessages(fect(
            Y ~ D,
            data   = d,
            index  = c("id", "time"),
            method = "fe",
            force  = "two-way",
            time.component.from = "nevertreated",
            se      = TRUE,
            vartype = "bootstrap",
            para.error = "wild",
            nboots  = 50,
            parallel = FALSE,
            CV = FALSE
        )))
    })
    ## vartype is bootstrap, para.error on fit object should be NULL
    ## (fect_boot para.error.resolved is only set in the parametric branch)
    expect_equal(fit$vartype, "bootstrap")
})

## -----------------------------------------------------------------------
## Test 9: para.error = "ar" works on fully-observed panel
## -----------------------------------------------------------------------

test_that("para.error = 'ar' runs on fully-observed panel", {
    skip_on_cran()
    d   <- make_panel()
    fit <- fect_para(d, para.error = "ar")
    expect_equal(fit$para.error, "ar")

    res <- estimand(fit, "att", "overall", window = c(1, 8), ci.method = "bca")
    expect_true(res$ci.lo < 3.0 && res$ci.hi > 3.0,
        info = sprintf("BCA CI [%.3f, %.3f] should contain 3.0", res$ci.lo, res$ci.hi))
})

## -----------------------------------------------------------------------
## Test 10: para.error = "auto" resolves to "ar" on missing-data panel
## -----------------------------------------------------------------------

test_that("para.error = 'auto' resolves to 'ar' on missing-data panel", {
    skip_on_cran()
    d   <- make_panel(add_missing = TRUE)
    ## na.rm = FALSE needed so missing cells are preserved in I
    fit <- suppressWarnings(suppressMessages(fect(
        Y ~ D,
        data   = d,
        index  = c("id", "time"),
        method = "fe",
        force  = "two-way",
        time.component.from = "nevertreated",
        se      = TRUE,
        vartype = "parametric",
        para.error = "auto",
        na.rm   = FALSE,
        nboots  = 100,
        parallel = FALSE,
        keep.sims = TRUE,
        CV = FALSE
    )))
    expect_equal(fit$para.error, "ar")
})
