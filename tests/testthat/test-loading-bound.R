## Integration tests for loading.bound argument in fect().
## See statsclaw-workspace/fect/runs/REQ-bounded-loadings/test-spec.md section I.

## ---------------------------------------------------------------------------
## Helper DGP: factor-structure panel with never-treated controls
## ---------------------------------------------------------------------------

make_factor_panel <- function(N = 100, TT = 30, Ntr = 30, tau = 3.0,
                              r = 2, seed = 42) {
    set.seed(seed)
    F_mat    <- matrix(rnorm(TT * r), TT, r)
    L_mat    <- matrix(rnorm(N * r), N, r)
    alpha_i  <- rnorm(N, 0, 1)
    xi_t     <- rnorm(TT, 0, 0.5)
    T0_vec   <- rep(Inf, N)
    if (Ntr > 0) {
        T0_vec[1:Ntr] <- sample(round(TT * 0.4):round(TT * 0.7), Ntr,
                                replace = TRUE)
    }
    Y_vec <- D_vec <- numeric(N * TT)
    id_vec <- time_vec <- integer(N * TT)
    idx <- 1
    for (i in seq_len(N)) {
        for (t in seq_len(TT)) {
            treated       <- (t >= T0_vec[i])
            D_vec[idx]    <- as.integer(treated)
            Y_vec[idx]    <- alpha_i[i] + xi_t[t] +
                             sum(F_mat[t, ] * L_mat[i, ]) +
                             tau * D_vec[idx] + rnorm(1, 0, 0.5)
            id_vec[idx]   <- i
            time_vec[idx] <- t
            idx           <- idx + 1
        }
    }
    data.frame(id = id_vec, time = time_vec, Y = Y_vec, D = D_vec)
}

common_fit <- function(df, ...) {
    suppressMessages(suppressWarnings(fect::fect(
        Y ~ D, data = df, index = c("id", "time"),
        method = "ife", r = 2, CV = FALSE, se = FALSE,
        parallel = FALSE, force = "two-way",
        time.component.from = "nevertreated", ...
    )))
}

## ---------------------------------------------------------------------------
## I.1: back-compat bit-identity with pinned dev@69bf243 output
## ---------------------------------------------------------------------------

test_that("I.1: default loading.bound = 'none' matches pinned dev@69bf243 output", {
    fixture_path <- test_path("fixtures", "bounded-loadings-backcompat-dev69bf243.rds")
    skip_if_not(file.exists(fixture_path), "Back-compat fixture not pinned.")
    pinned <- readRDS(fixture_path)

    df  <- make_factor_panel(N = 100, TT = 30, Ntr = 30, tau = 3.0, r = 2, seed = 42)
    set.seed(100)
    fit <- common_fit(df)
    expect_equal(fit$att.avg,   pinned$att.avg,   tolerance = 1e-10)
    expect_equal(fit$lambda.tr, pinned$lambda.tr, tolerance = 1e-10)
    expect_equal(fit$lambda.co, pinned$lambda.co, tolerance = 1e-10)
    expect_equal(fit$wgt.implied, pinned$wgt.implied, tolerance = 1e-10)
})

## ---------------------------------------------------------------------------
## I.2: simplex bound yields weights on simplex and lambda.tr = Lambda.co' %*% w
## ---------------------------------------------------------------------------

test_that("I.2: simplex bound produces lambda.tr in conv(lambda.co)", {
    df <- make_factor_panel(seed = 42)
    set.seed(100)
    fit <- common_fit(df, loading.bound = "simplex", gamma.loading = 1.0)

    expect_identical(fit$loading.bound, "simplex")
    expect_equal(fit$gamma.loading, 1.0, tolerance = 1e-12)

    for (i in seq_len(nrow(fit$lambda.tr))) {
        w_i <- fit$wgt.implied[i, ]
        expect_equal(sum(w_i), 1, tolerance = 1e-6)
        expect_true(all(w_i >= -1e-10))
        expect_equal(as.numeric(fit$lambda.tr[i, ]),
                     as.numeric(crossprod(fit$lambda.co, w_i)),
                     tolerance = 1e-6)
    }
})

## ---------------------------------------------------------------------------
## I.3: pointwise hull of fitted factor-component counterfactual
## ---------------------------------------------------------------------------

test_that("I.3: Y_hat(0) factor component lies in pointwise hull of F %*% lambda.co", {
    df <- make_factor_panel(seed = 42)
    set.seed(100)
    fit <- common_fit(df, loading.bound = "simplex", gamma.loading = 1.0)

    F_hat    <- fit$factor                       # TT x r
    Lam_co   <- fit$lambda.co                    # Nco x r
    ctrl_fit <- F_hat %*% t(Lam_co)              # TT x Nco
    Y_factor <- F_hat %*% t(fit$lambda.tr)       # TT x Ntr

    for (t in seq_len(nrow(ctrl_fit))) {
        lo <- min(ctrl_fit[t, ]); hi <- max(ctrl_fit[t, ])
        expect_true(all(Y_factor[t, ] >= lo - 1e-6 &
                        Y_factor[t, ] <= hi + 1e-6))
    }
})

## ---------------------------------------------------------------------------
## I.4: method restriction
## ---------------------------------------------------------------------------

test_that("I.4: loading.bound='simplex' errors with unsupported method/dispatch", {
    df <- make_factor_panel(N = 60, TT = 20, Ntr = 15, seed = 3)
    ## method = "mc" is independently blocked by the method check, but because
    ## the default `time.component.from = "notyettreated"`, either check may fire
    ## first. Both are acceptable errors indicating the configuration is invalid.
    expect_error(
        suppressMessages(suppressWarnings(fect(
            Y ~ D, data = df, index = c("id","time"),
            method = "mc", loading.bound = "simplex",
            CV = FALSE, se = FALSE, parallel = FALSE
        ))),
        regexp = "only supported for method|nevertreated"
    )
    ## method = "fe" is rewritten internally to method = "ife" with r = 0,
    ## so the effective failure mode is the time.component.from check.
    expect_error(
        suppressMessages(suppressWarnings(fect(
            Y ~ D, data = df, index = c("id","time"),
            method = "fe", loading.bound = "simplex",
            CV = FALSE, se = FALSE, parallel = FALSE
        ))),
        regexp = "only supported for method|nevertreated"
    )
})

test_that("I.4b: loading.bound='simplex' requires time.component.from='nevertreated'", {
    df <- make_factor_panel(seed = 5)
    expect_error(
        suppressMessages(suppressWarnings(fect(
            Y ~ D, data = df, index = c("id","time"),
            method = "ife", loading.bound = "simplex",
            time.component.from = "notyettreated",
            CV = FALSE, se = FALSE, parallel = FALSE, r = 2
        ))),
        regexp = "nevertreated"
    )
})

## ---------------------------------------------------------------------------
## I.5: argument validation
## ---------------------------------------------------------------------------

test_that("I.5: invalid loading.bound / gamma.loading / grid raise", {
    df <- make_factor_panel(seed = 7)
    expect_error(common_fit(df, loading.bound = "foo"),
                 regexp = "must be one of")
    expect_error(common_fit(df, loading.bound = "simplex", gamma.loading = -1),
                 regexp = "positive")
    expect_error(common_fit(df, loading.bound = "simplex",
                            gamma.loading.grid = c(1, -2, 3)),
                 regexp = "positive")
})

## ---------------------------------------------------------------------------
## I.6: gamma CV selects from grid
## ---------------------------------------------------------------------------

test_that("I.6: gamma.loading = NULL selects gamma from provided grid", {
    df <- make_factor_panel(seed = 9)
    grid <- c(0.1, 1.0, 10.0)
    set.seed(100)
    fit <- common_fit(df, loading.bound = "simplex",
                       gamma.loading = NULL, gamma.loading.grid = grid)
    expect_true(fit$gamma.loading %in% grid)
})

## ---------------------------------------------------------------------------
## I.7: reproducibility with fixed gamma
## ---------------------------------------------------------------------------

test_that("I.7: gamma fixed + same seed produces identical lambda.tr and wgt.implied", {
    df <- make_factor_panel(seed = 11)
    set.seed(100)
    fit1 <- common_fit(df, loading.bound = "simplex", gamma.loading = 2.5)
    set.seed(100)
    fit2 <- common_fit(df, loading.bound = "simplex", gamma.loading = 2.5)
    expect_equal(fit1$att.avg, fit2$att.avg, tolerance = 1e-10)
    expect_equal(fit1$lambda.tr, fit2$lambda.tr, tolerance = 1e-10)
    expect_equal(fit1$wgt.implied, fit2$wgt.implied, tolerance = 1e-10)
})

## ---------------------------------------------------------------------------
## I.8: diagnostic populated
## ---------------------------------------------------------------------------

test_that("I.8: loading.proj.resid is populated and sensible", {
    df <- make_factor_panel(seed = 13)
    set.seed(100)
    fit <- common_fit(df, loading.bound = "simplex", gamma.loading = 1.0)
    expect_false(is.null(fit$loading.proj.resid))
    expect_length(fit$loading.proj.resid, nrow(fit$lambda.tr))
    expect_true(all(fit$loading.proj.resid >= 0))
    expect_true(all(is.finite(fit$loading.proj.resid)))
})

## ---------------------------------------------------------------------------
## I.9: smoke test with bootstrap
## ---------------------------------------------------------------------------

test_that("I.9: fect() with se=TRUE runs with bounded loading (smoke)", {
    df <- make_factor_panel(N = 60, TT = 25, Ntr = 20, seed = 15)
    set.seed(100)
    fit <- suppressMessages(suppressWarnings(fect(
        Y ~ D, data = df, index = c("id","time"),
        method = "ife", r = 2, CV = FALSE, se = TRUE, nboots = 20,
        parallel = FALSE, force = "two-way",
        time.component.from = "nevertreated",
        loading.bound = "simplex", gamma.loading = 1.0
    )))
    expect_false(is.null(fit$est.att))
    expect_true(is.finite(fit$att.avg))
})

## ---------------------------------------------------------------------------
## I.10: loading.overlap plot type
## ---------------------------------------------------------------------------

test_that("I.10a: loading.overlap returns a ggplot for r >= 2", {
    df <- make_factor_panel(N = 80, TT = 25, Ntr = 20, r = 2, seed = 7)
    fit <- common_fit(df)
    p <- plot(fit, type = "loading.overlap")
    expect_s3_class(p, "ggplot")

    rd <- plot(fit, type = "loading.overlap", return.data = TRUE)
    expect_true(is.list(rd) || inherits(rd, "data.frame"))
})

test_that("I.10b: loading.overlap returns a ggplot for r == 1 (mirror histogram)", {
    df <- make_factor_panel(N = 80, TT = 25, Ntr = 20, r = 1, seed = 8)
    fit <- suppressMessages(suppressWarnings(fect::fect(
        Y ~ D, data = df, index = c("id", "time"),
        method = "ife", r = 1, CV = FALSE, se = FALSE,
        parallel = FALSE, force = "two-way",
        time.component.from = "nevertreated"
    )))
    p <- plot(fit, type = "loading.overlap")
    expect_s3_class(p, "ggplot")
})

test_that("I.10c: loading.overlap errors when r == 0", {
    df  <- make_factor_panel(N = 60, TT = 20, Ntr = 15, r = 1, seed = 9)
    fit <- suppressMessages(suppressWarnings(fect::fect(
        Y ~ D, data = df, index = c("id", "time"),
        method = "fe", r = 0, CV = FALSE, se = FALSE,
        parallel = FALSE, force = "two-way"
    )))
    expect_error(plot(fit, type = "loading.overlap"), regexp = "r >= 1")
})
