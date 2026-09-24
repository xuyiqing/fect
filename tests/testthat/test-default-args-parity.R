## fect() is an S3 generic. A formula call goes to fect.formula(), which passes
## every argument on to fect.default(). A call that names the variables as
## strings (Y = "y", D = "d", ...) goes straight to fect.default(), and
## UseMethod() does not carry the generic's defaults along. So fect.default()
## must declare the same defaults itself. Before 2.4.6 it did not: cv.method was
## "all_units" (block CV) instead of "rolling" (the default since v2.3.0), and
## nlambda was 0 instead of 10, so method = "mc" without lambda stopped with
## '"nlambda" option misspecified.'

skip_on_cran()

## Deparsed default of every argument, named by argument ("" = no default).
.fect_defaults <- function(f) {
    fm <- formals(f)
    vapply(names(fm), function(nm) paste(deparse(fm[[nm]]), collapse = " "),
           character(1))
}

## Two factors (the second one weaker), one covariate, staggered adoption.
## With seed = 10, block CV selects r = 1 and rolling CV selects r = 2 (as of
## 2.4.6), so a call that silently fell back to block CV changes r.cv.
.parity_panel <- function(seed = 10, N = 30, TT = 16, Ntr = 10, T0 = 10) {
    set.seed(seed)
    alpha <- rnorm(N)
    xi <- rnorm(TT)
    lam <- cbind(rnorm(N), 0.5 * rnorm(N))
    fac <- matrix(rnorm(TT * 2), TT, 2)
    X <- matrix(rnorm(N * TT), TT, N)
    D <- matrix(0, TT, N)
    for (j in seq_len(Ntr)) D[(T0 + 1 + (j %% 3)):TT, j] <- 1
    Y <- outer(xi, rep(1, N)) + outer(rep(1, TT), alpha) + fac %*% t(lam) +
        0.5 * X + 3 * D + matrix(rnorm(N * TT), TT, N)
    data.frame(id = rep(seq_len(N), each = TT), time = rep(seq_len(TT), N),
               Y = c(Y), D = c(D), X = c(X))
}

test_that("fect(), fect.formula() and fect.default() declare the same defaults", {
    generic <- .fect_defaults(fect::fect)
    expect_identical(.fect_defaults(fect:::fect.formula), generic)
    expect_identical(.fect_defaults(fect:::fect.default), generic)
})

test_that("formula and non-formula calls use the same CV method and select the same r", {
    d <- .parity_panel()
    ## Seed the session right before each call so both calls draw the same CV
    ## folds. (Passing seed = 1 does the same since 2.4.6; set.seed() keeps
    ## this test about calling styles only.)
    set.seed(1)
    out_formula <- suppressMessages(fect(
        Y ~ D + X, data = d, index = c("id", "time"),
        method = "ife", CV = TRUE, r = c(0, 3), parallel = FALSE))
    set.seed(1)
    msgs <- capture_messages(out_default <- fect(
        Y = "Y", D = "D", X = "X", data = d, index = c("id", "time"),
        method = "ife", CV = TRUE, r = c(0, 3), parallel = FALSE))

    ## Same folds and the same scores: both calls ran the same CV method.
    expect_equal(out_default$CV.out.ife, out_formula$CV.out.ife)
    expect_identical(out_default$r.cv, out_formula$r.cv)
    ## The old default also printed a deprecation note for "all_units", an
    ## option the user never set.
    expect_false(any(grepl("all_units", msgs, fixed = TRUE)))
})

test_that("a non-formula call with method = \"mc\" and no lambda cross-validates lambda", {
    d <- .parity_panel()
    set.seed(1)
    out_formula <- suppressMessages(fect(
        Y ~ D + X, data = d, index = c("id", "time"),
        method = "mc", parallel = FALSE))
    set.seed(1)
    out_default <- suppressMessages(fect(
        Y = "Y", D = "D", X = "X", data = d, index = c("id", "time"),
        method = "mc", parallel = FALSE))

    expect_equal(out_default$CV.out.mc, out_formula$CV.out.mc)
    expect_equal(out_default$lambda.cv, out_formula$lambda.cv)
})
