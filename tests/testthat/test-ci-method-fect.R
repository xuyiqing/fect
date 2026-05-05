## Tests for the v2.4.2 ci.method = c("normal", "basic") arg on fect(),
## and the soft-deprecated quantile.CI legacy arg.
##
## What we're verifying:
##   1. Default behavior is byte-equal to the pre-v2.4.2 default
##      (quantile.CI = FALSE -> normal Wald CI).
##   2. ci.method = "normal" is byte-equal to legacy quantile.CI = FALSE.
##   3. ci.method = "basic" is byte-equal to legacy quantile.CI = TRUE.
##   4. Supplying quantile.CI explicitly emits a deprecation warning.
##   5. NOT supplying quantile.CI (NULL sentinel) emits no warning under the
##      modern API.
##   6. ci.method = "basic" with nboots < 1000 emits a tail-CI replicate
##      warning at fit time (mirrors the estimand .check_tail_ci_replicates
##      gate).
##   7. ci.method = "bca" / "bc" / "percentile" hard-error with a message
##      pointing to estimand() for the full 5-method surface.
##   8. ci.method = "garbage" hard-errors with a clear message.

library(testthat)

data(simdata, package = "fect")

base_call <- function(nboots = 50, parallel = FALSE, ...) {
    fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
         method = "fe", force = "two-way", se = TRUE, nboots = nboots,
         parallel = parallel, seed = 42L, ...)
}

# ---- 1. Default behavior is byte-equal to legacy quantile.CI = FALSE -------

test_that("default ci.method yields the legacy quantile.CI = FALSE result", {
    out_default <- base_call()
    out_legacy_false <- suppressWarnings(base_call(quantile.CI = FALSE))
    expect_equal(out_default$est.avg, out_legacy_false$est.avg)
    expect_equal(out_default$est.att, out_legacy_false$est.att)
})

# ---- 2. ci.method = "normal" byte-equal to legacy FALSE --------------------

test_that("ci.method = 'normal' is byte-equal to quantile.CI = FALSE", {
    out_normal <- base_call(ci.method = "normal")
    out_legacy <- suppressWarnings(base_call(quantile.CI = FALSE))
    expect_equal(out_normal$est.avg, out_legacy$est.avg)
    expect_equal(out_normal$est.att, out_legacy$est.att)
})

# ---- 3. ci.method = "basic" byte-equal to legacy TRUE ----------------------

test_that("ci.method = 'basic' is byte-equal to quantile.CI = TRUE", {
    out_basic <- suppressWarnings(base_call(ci.method = "basic", nboots = 1000))
    out_legacy <- suppressWarnings(base_call(quantile.CI = TRUE, nboots = 1000))
    expect_equal(out_basic$est.avg, out_legacy$est.avg)
    expect_equal(out_basic$est.att, out_legacy$est.att)
})

# ---- 4. Supplying quantile.CI explicitly fires deprecation warning ---------

test_that("quantile.CI = TRUE fires a deprecation warning that mentions ci.method", {
    expect_warning(
        base_call(quantile.CI = TRUE, nboots = 1000),
        regexp = "quantile\\.CI.*deprecated.*ci\\.method"
    )
})

test_that("quantile.CI = FALSE also fires the deprecation warning", {
    ## "Soft-deprecate" means any user-supplied use of the legacy arg warns,
    ## even when supplying the legacy default value.
    expect_warning(
        base_call(quantile.CI = FALSE),
        regexp = "quantile\\.CI.*deprecated"
    )
})

# ---- 5. NULL sentinel: no warning under modern API -------------------------

test_that("modern API (no quantile.CI supplied) emits no quantile.CI warning", {
    expect_silent_about <- function(expr, regexp) {
        warned <- FALSE
        withCallingHandlers(
            expr,
            warning = function(w) {
                if (grepl(regexp, conditionMessage(w))) warned <<- TRUE
                invokeRestart("muffleWarning")
            }
        )
        expect_false(warned)
    }
    expect_silent_about(base_call(), "quantile\\.CI")
    expect_silent_about(base_call(ci.method = "normal"), "quantile\\.CI")
})

# ---- 6. nboots < 1000 + ci.method = "basic" warns at fit time --------------

test_that("ci.method = 'basic' with nboots < 1000 fires the tail-CI warning", {
    expect_warning(
        base_call(ci.method = "basic", nboots = 50),
        regexp = "tail quantiles.*1000|Efron 1987"
    )
})

expect_silent_about <- function(expr, regexp) {
    warned <- FALSE
    withCallingHandlers(
        expr,
        warning = function(w) {
            if (grepl(regexp, conditionMessage(w))) warned <<- TRUE
            invokeRestart("muffleWarning")
        }
    )
    expect_false(warned)
}

test_that("ci.method = 'basic' with nboots >= 1000 emits no tail-CI warning", {
    expect_silent_about(
        base_call(ci.method = "basic", nboots = 1000),
        "tail quantiles|Efron 1987"
    )
})

test_that("ci.method = 'normal' never fires the tail-CI warning", {
    expect_silent_about(
        base_call(ci.method = "normal", nboots = 50),
        "tail quantiles"
    )
})

# ---- 7. ci.method in {bca, bc, percentile} hard-errors ---------------------

test_that("ci.method = 'bca' errors with a message pointing to estimand()", {
    expect_error(
        base_call(ci.method = "bca"),
        regexp = "bca.*not supported in fect.*estimand"
    )
})

test_that("ci.method = 'bc' errors with a message pointing to estimand()", {
    expect_error(
        base_call(ci.method = "bc"),
        regexp = "bc.*not supported in fect.*estimand"
    )
})

test_that("ci.method = 'percentile' errors with a message pointing to estimand()", {
    expect_error(
        base_call(ci.method = "percentile"),
        regexp = "percentile.*not supported in fect.*estimand"
    )
})

# ---- 8. ci.method validation -----------------------------------------------

test_that("ci.method = 'garbage' fails fast", {
    expect_error(
        base_call(ci.method = "garbage"),
        regexp = "ci\\.method.*must be one of"
    )
})

# ---- 9. ci.method = "basic" + vartype = "jackknife" -> error ---------------

test_that("ci.method = 'basic' + vartype = 'jackknife' hard-errors", {
    expect_error(
        suppressWarnings(base_call(vartype = "jackknife", ci.method = "basic")),
        regexp = "ci\\.method = \"basic\".*not supported.*jackknife"
    )
})

# ---- 10. ci.method = "basic" + parametric -> location-shift fix ------------
#
# fect's built-in CI machinery applies a location-shift fix at the avg-level
# and per-event-time CI sites so that ci.method = "basic" works on parametric
# fits.  The fect-level result should match estimand(fit, "att", "basic")
# byte-equally on the avg-level CI.

test_that("ci.method = 'basic' + parametric matches estimand at avg level", {
    skip_if_not_installed("fect")
    # Build a simple no-reversal panel (parametric requires no reversal)
    set.seed(1001L); N <- 30; TT <- 20; T0 <- 12; Ntr <- 5
    alpha_i <- rnorm(N); xi_t <- rnorm(TT)
    e <- matrix(rnorm(TT * N), TT, N)
    D <- matrix(0L, TT, N); D[(T0 + 1):TT, 1:Ntr] <- 1L
    Y <- outer(xi_t, rep(1, N)) + outer(rep(1, TT), alpha_i) + 3 * D + e
    df <- data.frame(id = rep(1:N, each = TT), time = rep(1:TT, N),
                     Y = c(Y), D = c(D))
    fit <- suppressWarnings(fect(
        Y ~ D, data = df, index = c("id", "time"),
        method = "fe", force = "two-way", CV = FALSE,
        se = TRUE, vartype = "parametric", para.error = "auto",
        ci.method = "basic",
        nboots = 200, parallel = FALSE,
        time.component.from = "nevertreated",
        keep.sims = TRUE, seed = 42L
    ))
    e <- suppressWarnings(estimand(fit, type = "att", by = "overall",
                                    ci.method = "basic"))
    expect_equal(unname(fit$est.avg[1, "CI.lower"]), e$ci.lo, tolerance = 1e-10)
    expect_equal(unname(fit$est.avg[1, "CI.upper"]), e$ci.hi, tolerance = 1e-10)
})

# ---- 9. vartype = "jackknife" + cl emits an "ignored" warning --------------
#
# fect's jackknife is leave-one-unit-out and does not consult cl.  Combining
# vartype = "jackknife" with a non-NULL cl silently produced unit-level SEs.
# The user-facing warning routes such callers to vartype = "bootstrap" + cl
# for cluster-aware inference.

test_that("vartype = 'jackknife' with cl warns that cl is ignored", {
    d <- simdata
    d$cl <- (d$id - 1) %/% 2 + 1   # 2 units per cluster
    expect_warning(
        suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"),
                               method = "fe", force = "two-way",
                               vartype = "jackknife", cl = "cl",
                               se = TRUE, parallel = FALSE, seed = 42L)),
        "cl argument is ignored"
    )
})

test_that("vartype = 'jackknife' without cl does not emit the cl warning", {
    seen <- character(0)
    withCallingHandlers(
        suppressMessages(fect(Y ~ D, data = simdata, index = c("id", "time"),
                               method = "fe", force = "two-way",
                               vartype = "jackknife",
                               se = TRUE, parallel = FALSE, seed = 42L)),
        warning = function(w) {
            seen <<- c(seen, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    expect_false(any(grepl("cl argument is ignored", seen, fixed = TRUE)))
})

test_that("vartype = 'bootstrap' with cl does not emit the jackknife warning", {
    d <- simdata
    d$cl <- (d$id - 1) %/% 2 + 1
    seen <- character(0)
    withCallingHandlers(
        suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"),
                               method = "fe", force = "two-way",
                               vartype = "bootstrap", cl = "cl",
                               nboots = 50, se = TRUE,
                               parallel = FALSE, seed = 42L)),
        warning = function(w) {
            seen <<- c(seen, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    expect_false(any(grepl("cl argument is ignored", seen, fixed = TRUE)))
})
