## ---------------------------------------------------------
## Tests derived from Quarto book documentation claims
## Each test validates a specific assertion made in the user-facing docs.
## ---------------------------------------------------------

## =========================================================
## Shared fixtures (computed once, reused across test blocks)
## =========================================================

## DGP helper: two-way FE + optional factors + optional covariates
make_panel <- function(N = 40, TT = 20, T0 = 12, Ntr = 12,
                       tau = 3.0, r = 0, beta_x = 0,
                       mean_shift = 0, seed = 8001,
                       stagger = FALSE, reversals = FALSE,
                       add_Z = FALSE, add_group = FALSE) {
    set.seed(seed)
    alpha_i <- rnorm(N, 0, 2)
    xi_t <- rnorm(TT, 0, 1)
    ## factors
    FL <- matrix(0, TT, N)
    if (r > 0) {
        Fmat <- matrix(rnorm(TT * r), TT, r)
        Lmat <- matrix(rnorm(N * r), N, r)
        FL <- Fmat %*% t(Lmat)
    }
    ## treatment assignment
    D <- matrix(0L, TT, N)
    tr_units <- 1:Ntr
    if (stagger) {
        ## 3 cohorts
        for (i in tr_units) {
            cohort <- ((i - 1) %% 3) + 1
            t0_i <- T0 + (cohort - 1) * 2
            if (t0_i < TT) D[(t0_i + 1):TT, i] <- 1L
        }
    } else {
        D[(T0 + 1):TT, tr_units] <- 1L
    }
    if (reversals) {
        ## flip last 2 periods off for first 3 treated
        for (i in 1:min(3, Ntr)) {
            D[(TT - 1):TT, i] <- 0L
        }
    }
    ## covariates
    X1 <- matrix(rnorm(N * TT, 0, 1), TT, N)
    eps <- matrix(rnorm(N * TT, 0, 1), TT, N)
    Y <- mean_shift +
        outer(xi_t, rep(1, N)) +
        outer(rep(1, TT), alpha_i) +
        FL +
        beta_x * X1 +
        tau * D +
        eps
    df <- data.frame(
        id = rep(1:N, each = TT),
        time = rep(1:TT, N),
        Y = as.vector(Y),
        D = as.vector(D),
        X1 = as.vector(X1)
    )
    if (add_Z) {
        ## time-invariant covariate (same across all periods for a unit)
        Z_vals <- rnorm(N, 0, 1)
        df$Z <- rep(Z_vals, each = TT)
    }
    if (add_group) {
        ## group variable: 3 groups by unit
        df$grp <- rep(((1:N - 1) %% 3) + 1, each = TT)
    }
    df
}


## =========================================================
## Section A: Method Equivalences (Book: index.qmd, 04-gsynth.Rmd)
## =========================================================

test_that("A1: gsynth == ife + time.component.from='nevertreated' (exact)", {
    skip_on_cran()
    d <- make_panel(N = 50, TT = 20, T0 = 12, Ntr = 10, r = 2, seed = 8010)
    out_g <- fect(Y ~ D + X1, data = d, index = c("id", "time"),
                  method = "gsynth", r = 2, se = FALSE, CV = FALSE)
    out_i <- fect(Y ~ D + X1, data = d, index = c("id", "time"),
                  method = "ife", time.component.from = "nevertreated",
                  r = 2, se = FALSE, CV = FALSE)
    expect_equal(out_g$att.avg, out_i$att.avg, tolerance = 1e-8)
    expect_equal(as.vector(out_g$eff), as.vector(out_i$eff), tolerance = 1e-8)
})

test_that("A2: IFE(r=0) == FE (exact equivalence)", {
    skip_on_cran()
    d <- make_panel(N = 30, TT = 15, T0 = 10, Ntr = 10, r = 0, seed = 8020)
    out_fe <- fect(Y ~ D + X1, data = d, index = c("id", "time"),
                   method = "fe", se = FALSE, CV = FALSE)
    out_ife <- fect(Y ~ D + X1, data = d, index = c("id", "time"),
                    method = "ife", r = 0, se = FALSE, CV = FALSE)
    expect_equal(out_fe$att.avg, out_ife$att.avg, tolerance = 1e-6)
    expect_equal(as.vector(out_fe$eff), as.vector(out_ife$eff), tolerance = 1e-6)
})


## =========================================================
## Section B: force Parameter Variations (Book: 02-fect.Rmd)
## =========================================================

test_that("B1: force='unit' runs without error", {
    skip_on_cran()
    d <- make_panel(N = 30, TT = 15, T0 = 10, Ntr = 10, seed = 8030)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", force = "unit", se = FALSE, CV = FALSE)
    expect_true(is.numeric(out$att.avg))
    expect_false(is.na(out$att.avg))
})

test_that("B2: force='time' runs without error", {
    skip_on_cran()
    d <- make_panel(N = 30, TT = 15, T0 = 10, Ntr = 10, seed = 8031)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", force = "time", se = FALSE, CV = FALSE)
    expect_true(is.numeric(out$att.avg))
    expect_false(is.na(out$att.avg))
})

test_that("B3: force='none' runs without error", {
    skip_on_cran()
    d <- make_panel(N = 30, TT = 15, T0 = 10, Ntr = 10, seed = 8032)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", force = "none", se = FALSE, CV = FALSE)
    expect_true(is.numeric(out$att.avg))
    expect_false(is.na(out$att.avg))
})

test_that("B4: force variations work for IFE with r > 0", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12, r = 1, seed = 8033)
    for (f in c("unit", "time", "none", "two-way")) {
        out <- fect(Y ~ D, data = d, index = c("id", "time"),
                    method = "ife", r = 1, force = f, se = FALSE, CV = FALSE)
        expect_true(is.numeric(out$att.avg),
                    info = paste("force =", f))
    }
})


## =========================================================
## Section C: Diagnostic Tests (Book: 02-fect.Rmd)
## =========================================================

test_that("C1: placeboTest produces p-value and test output", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12, seed = 8040)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", se = TRUE, nboots = 50,
                placeboTest = TRUE, placebo.period = c(-2, 0),
                CV = FALSE)
    ## Book: placebo p-value should exist
    expect_true(!is.null(out$placebo.p))
    expect_true(is.numeric(out$placebo.p))
})

test_that("C2: carryoverTest produces test output", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12,
                    reversals = TRUE, seed = 8041)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", se = TRUE, nboots = 50,
                carryoverTest = TRUE, carryover.period = c(1, 2),
                CV = FALSE)
    expect_true(!is.null(out$carryover.p))
    expect_true(is.numeric(out$carryover.p))
})

test_that("C3: loo=TRUE produces pre-trend leave-one-out estimates", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12, seed = 8042)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", se = TRUE, nboots = 50,
                loo = TRUE, CV = FALSE)
    ## LOO test output should exist
    expect_true(!is.null(out$loo.test.out))
    expect_true(!is.null(out$loo.test.out$f.stat))
    expect_true(!is.null(out$loo.test.out$f.p))
})

test_that("C4: equivalence test output (pre-trend F-test)", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12, seed = 8043)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", se = TRUE, nboots = 50, loo = TRUE, CV = FALSE)
    ## F-test lives inside loo.test.out
    expect_true(!is.null(out$loo.test.out$f.p))
    expect_true(is.numeric(out$loo.test.out$f.p))
    ## Book: p > 0.05 means good pre-trend fitting; for clean DGP expect pass
    expect_true(out$loo.test.out$f.p > 0.01,
                info = "Pre-trend F-test should not reject for clean DGP")
})


## =========================================================
## Section D: Output Structure (Book: 02-fect.Rmd, 04-gsynth.Rmd)
## =========================================================

test_that("D1: est.att structure (periods × 6 columns)", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12, seed = 8050)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", se = TRUE, nboots = 50, CV = FALSE)
    expect_true(!is.null(out$est.att))
    expect_true(is.matrix(out$est.att) || is.data.frame(out$est.att))
    ## Book: columns include ATT, S.E., CI.lower, CI.upper, p.value, count
    cnames <- colnames(out$est.att)
    expect_true("ATT" %in% cnames)
    expect_true("S.E." %in% cnames)
    expect_true("CI.lower" %in% cnames)
    expect_true("CI.upper" %in% cnames)
    expect_true("count" %in% cnames)
})

test_that("D2: est.avg vs est.avg.unit differ with unequal exposure", {
    skip_on_cran()
    ## staggered adoption → unequal treated-period counts
    d <- make_panel(N = 40, TT = 20, T0 = 10, Ntr = 12,
                    stagger = TRUE, seed = 8051)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", se = FALSE, CV = FALSE)
    ## Book: est.avg weights observations equally; est.avg.unit weights units equally
    ## With staggered adoption these should differ
    expect_true(is.numeric(out$att.avg))
    expect_true(!is.null(out$att.avg.unit))
    ## They won't be exactly equal due to stagger
    expect_false(isTRUE(all.equal(out$att.avg, out$att.avg.unit)),
                 info = "est.avg and est.avg.unit should differ under staggered adoption")
})

test_that("D3: eff.boot dimensions = periods × treated × nboots", {
    skip_on_cran()
    d <- make_panel(N = 30, TT = 15, T0 = 10, Ntr = 8, seed = 8052)
    nboots <- 30
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", se = TRUE, nboots = nboots,
                keep.sims = TRUE, CV = FALSE)
    ## Book: eff.boot is array(periods × treated × nboots)
    expect_true(!is.null(out$eff.boot))
    expect_true(is.array(out$eff.boot))
    dims <- dim(out$eff.boot)
    expect_equal(length(dims), 3)
    ## dim 3 should be nboots
    expect_equal(dims[3], nboots)
})

test_that("D4: est.beta present with covariates", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12,
                    beta_x = 1.5, seed = 8053)
    out <- fect(Y ~ D + X1, data = d, index = c("id", "time"),
                method = "fe", se = FALSE, CV = FALSE)
    ## Book: beta reports time-varying covariate coefficients
    expect_true(!is.null(out$beta))
    expect_true(is.numeric(out$beta))
})

test_that("D5: sigma2 is estimated error variance", {
    skip_on_cran()
    d <- make_panel(N = 50, TT = 20, T0 = 12, Ntr = 15, seed = 8054)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", se = FALSE, CV = FALSE)
    ## DGP has eps ~ N(0,1) so sigma2 should be close to 1
    expect_true(!is.null(out$sigma2))
    expect_true(out$sigma2 > 0)
    expect_true(out$sigma2 < 3,
                info = "sigma2 should be in reasonable range for sd=1 DGP")
})

test_that("D6: wgt.implied dimensions = Nco × Ntr for gsynth", {
    skip_on_cran()
    d <- make_panel(N = 50, TT = 20, T0 = 12, Ntr = 10, r = 1, seed = 8055)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "gsynth", r = 1, se = FALSE, CV = FALSE)
    if (!is.null(out$wgt.implied)) {
        ## Book: dimensions = N_co × N_tr
        expect_equal(ncol(out$wgt.implied), 10)  ## Ntr
        expect_equal(nrow(out$wgt.implied), 40)   ## N - Ntr
        ## Book: weights can be positive and negative
        expect_true(any(out$wgt.implied > 0))
    }
})


## =========================================================
## Section E: Edge Cases from Book (02-fect.Rmd, 04-gsynth.Rmd)
## =========================================================

test_that("E1: single treated unit + parametric bootstrap", {
    skip_on_cran()
    d <- make_panel(N = 30, TT = 20, T0 = 12, Ntr = 1, r = 1, seed = 8060)
    ## Book: parametric bootstrap works even with Ntr = 1
    out <- suppressWarnings(fect(
        Y ~ D, data = d, index = c("id", "time"),
        method = "gsynth", r = 1, se = TRUE,
        vartype = "parametric", nboots = 30, CV = FALSE
    ))
    expect_true(is.numeric(out$att.avg))
    expect_false(is.na(out$att.avg))
    ## CI should exist
    expect_true(!is.null(out$est.att))
})

test_that("E2: missing Y allowed, missing D errors", {
    skip_on_cran()
    d <- make_panel(N = 30, TT = 15, T0 = 10, Ntr = 8, seed = 8061)
    ## Introduce missing Y
    d$Y[sample(nrow(d), 10)] <- NA
    out <- suppressWarnings(fect(
        Y ~ D, data = d, index = c("id", "time"),
        method = "fe", se = FALSE, CV = FALSE, na.rm = TRUE
    ))
    expect_true(is.numeric(out$att.avg))

    ## Missing D should error
    d2 <- make_panel(N = 30, TT = 15, T0 = 10, Ntr = 8, seed = 8062)
    d2$D[5] <- NA
    expect_error(
        fect(Y ~ D, data = d2, index = c("id", "time"),
             method = "fe", se = FALSE, CV = FALSE),
        info = "Missing D should produce an error"
    )
})

test_that("E3: unbalanced panel (missing time periods) works", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12, seed = 8063)
    ## Drop ~10% of observations to create unbalanced panel
    keep <- sample(nrow(d), round(nrow(d) * 0.9))
    d_unbal <- d[keep, ]
    out <- suppressWarnings(fect(
        Y ~ D, data = d_unbal, index = c("id", "time"),
        method = "ife", r = 1, se = FALSE, CV = FALSE
    ))
    expect_true(is.numeric(out$att.avg))
    expect_false(is.na(out$att.avg))
})

test_that("E4: balance.period restricts included units", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 10, Ntr = 12,
                    stagger = TRUE, seed = 8064)
    out_full <- fect(Y ~ D, data = d, index = c("id", "time"),
                     method = "fe", se = FALSE, CV = FALSE)
    out_bal <- fect(Y ~ D, data = d, index = c("id", "time"),
                    method = "fe", se = FALSE, CV = FALSE,
                    balance.period = c(-3, 3))
    ## Balanced subset should have fewer or equal treated obs
    expect_true(sum(out_bal$count, na.rm = TRUE) <=
                sum(out_full$count, na.rm = TRUE))
})

test_that("E5: min.T0 drops units with few control observations", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12, seed = 8065)
    out5 <- fect(Y ~ D, data = d, index = c("id", "time"),
                 method = "fe", se = FALSE, CV = FALSE, min.T0 = 5)
    out10 <- fect(Y ~ D, data = d, index = c("id", "time"),
                  method = "fe", se = FALSE, CV = FALSE, min.T0 = 10)
    ## Stricter min.T0 should drop more or equal treated units
    expect_true(sum(out10$count, na.rm = TRUE) <=
                sum(out5$count, na.rm = TRUE))
})

test_that("E6: seed produces reproducible results", {
    skip_on_cran()
    d <- make_panel(N = 30, TT = 15, T0 = 10, Ntr = 8, seed = 8066)
    out1 <- fect(Y ~ D, data = d, index = c("id", "time"),
                 method = "fe", se = TRUE, nboots = 30, seed = 1234, CV = FALSE)
    out2 <- fect(Y ~ D, data = d, index = c("id", "time"),
                 method = "fe", se = TRUE, nboots = 30, seed = 1234, CV = FALSE)
    expect_equal(out1$att.avg, out2$att.avg, tolerance = 1e-10)
    expect_equal(out1$est.att, out2$est.att, tolerance = 1e-10)
})


## =========================================================
## Section F: MC Method (Book: 04-gsynth.Rmd cheatsheet)
## =========================================================

test_that("F1: method='mc' basic functionality", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12, r = 1, seed = 8070)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "mc", se = FALSE, CV = FALSE)
    expect_true(is.numeric(out$att.avg))
    expect_false(is.na(out$att.avg))
    ## ATT should be in reasonable range of true tau=3.0
    expect_true(abs(out$att.avg - 3.0) < 2.0)
})

test_that("F2: mc rejects time.component.from='nevertreated'", {
    skip_on_cran()
    d <- make_panel(N = 30, TT = 15, T0 = 10, Ntr = 8, seed = 8071)
    expect_error(
        fect(Y ~ D, data = d, index = c("id", "time"),
             method = "mc", time.component.from = "nevertreated",
             se = FALSE, CV = FALSE)
    )
})

test_that("F3: mc with CV selects lambda", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12, r = 1, seed = 8072)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "mc", se = FALSE, CV = TRUE)
    ## Should have selected a lambda
    expect_true(!is.null(out$lambda.cv) || !is.null(out$lambda))
})


## =========================================================
## Section G: CFE Features (Book: 07-cfe.Rmd)
## =========================================================

test_that("G1: CFE with unit-specific linear time trend (Q.type)", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12, seed = 8080)
    out <- fect(Y ~ D + X1, data = d, index = c("id", "time"),
                method = "cfe", r = 0, se = FALSE, CV = FALSE,
                Q.type = "linear")
    expect_true(is.numeric(out$att.avg))
    expect_false(is.na(out$att.avg))
})

test_that("G2: CFE with time-invariant covariate Z and gamma", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12,
                    add_Z = TRUE, seed = 8081)
    ## Create time group for gamma: early vs late
    d$time_group <- ifelse(d$time <= 10, 1, 2)
    out <- fect(Y ~ D + X1, data = d, index = c("id", "time"),
                method = "cfe", r = 0, se = FALSE, CV = FALSE,
                Z = "Z", gamma = "time_group")
    expect_true(is.numeric(out$att.avg))
    expect_false(is.na(out$att.avg))
})

test_that("G3: CFE with extra group FE via index[3]", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12,
                    add_group = TRUE, seed = 8082)
    ## Book: extra FE specified via index with 3+ elements
    out <- fect(Y ~ D + X1, data = d, index = c("id", "time", "grp"),
                method = "cfe", r = 0, se = FALSE, CV = FALSE)
    expect_true(is.numeric(out$att.avg))
    expect_false(is.na(out$att.avg))
})

test_that("G4: CFE defaults to time.component.from='notyettreated'", {
    skip_on_cran()
    d <- make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12, seed = 8083)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "cfe", r = 0, se = FALSE, CV = FALSE)
    ## Book: CFE defaults to notyettreated, unlike gsynth
    expect_equal(out$time.component.from, "notyettreated")
})

test_that("G5: CFE with r > 0 adds latent factors", {
    skip_on_cran()
    d <- make_panel(N = 50, TT = 20, T0 = 12, Ntr = 10,
                    r = 2, seed = 8084)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "cfe", r = 2, se = FALSE, CV = FALSE,
                time.component.from = "nevertreated")
    expect_true(is.numeric(out$att.avg))
    ## ATT should be in reasonable range
    expect_true(abs(out$att.avg - 3.0) < 2.0)
})


## =========================================================
## Section H: Inference Methods (Book: cheatsheet, 02-fect.Rmd)
## =========================================================

test_that("H1: vartype='parametric' only for gsynth/nevertreated", {
    skip_on_cran()
    d <- make_panel(N = 50, TT = 20, T0 = 12, Ntr = 10, r = 1, seed = 8090)
    ## Book: parametric specific to gsynth
    out <- suppressWarnings(fect(Y ~ D, data = d, index = c("id", "time"),
                method = "gsynth", r = 1, se = TRUE,
                vartype = "parametric", nboots = 30, CV = FALSE))
    expect_true(!is.null(out$est.att))
    ## Also works for ife with time.component.from='nevertreated' (Phase 3 unlocked)
    out2 <- suppressWarnings(fect(Y ~ D, data = d, index = c("id", "time"),
                 method = "ife", time.component.from = "nevertreated",
                 r = 1, se = TRUE, vartype = "parametric",
                 nboots = 30, CV = FALSE))
    expect_true(!is.null(out2$est.att))
})

test_that("H2: vartype='jackknife' works for FE and IFE", {
    skip_on_cran()
    d <- make_panel(N = 30, TT = 15, T0 = 10, Ntr = 8, seed = 8091)
    out_fe <- fect(Y ~ D, data = d, index = c("id", "time"),
                   method = "fe", se = TRUE, vartype = "jackknife", CV = FALSE)
    expect_true(!is.null(out_fe$est.att))

    out_ife <- fect(Y ~ D, data = d, index = c("id", "time"),
                    method = "ife", r = 0, se = TRUE,
                    vartype = "jackknife", CV = FALSE)
    expect_true(!is.null(out_ife$est.att))
})

test_that("H3: nboots default is 200", {
    skip_on_cran()
    ## Check default via formals
    defs <- formals(fect)
    expect_equal(defs$nboots, 200)
})


## =========================================================
## Section I: Plot Types (Book: 03-plots.Rmd, 04-gsynth.Rmd)
## =========================================================

test_that("I1: gap plot (default) returns ggplot", {
    skip_on_cran()
    d <- make_panel(N = 30, TT = 15, T0 = 10, Ntr = 8, seed = 8100)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", se = TRUE, nboots = 30, CV = FALSE)
    p <- plot(out, type = "gap")
    expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("I2: counterfactual plot returns ggplot", {
    skip_on_cran()
    d <- make_panel(N = 30, TT = 15, T0 = 10, Ntr = 8, seed = 8101)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", se = TRUE, nboots = 30, CV = FALSE)
    p <- plot(out, type = "ct", id = 1)
    expect_true(inherits(p, "gg") || inherits(p, "ggplot") || is.list(p))
})

test_that("I3: equivalence plot returns ggplot", {
    skip_on_cran()
    d <- make_panel(N = 30, TT = 15, T0 = 10, Ntr = 8, seed = 8102)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", se = TRUE, nboots = 30, CV = FALSE)
    p <- plot(out, type = "equiv")
    expect_true(inherits(p, "gg") || inherits(p, "ggplot") || is.list(p))
})

test_that("I4: status plot returns ggplot", {
    skip_on_cran()
    d <- make_panel(N = 30, TT = 15, T0 = 10, Ntr = 8, seed = 8103)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "fe", se = FALSE, CV = FALSE)
    p <- plot(out, type = "status")
    expect_true(inherits(p, "gg") || inherits(p, "ggplot") || is.list(p))
})

test_that("I5: gsynth factor and loading plots", {
    skip_on_cran()
    d <- make_panel(N = 50, TT = 20, T0 = 12, Ntr = 10, r = 2, seed = 8104)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "gsynth", r = 2, se = FALSE, CV = FALSE)
    pf <- plot(out, type = "factors")
    expect_true(inherits(pf, "gg") || inherits(pf, "ggplot") || is.list(pf))
    pl <- plot(out, type = "loadings")
    expect_true(inherits(pl, "gg") || inherits(pl, "ggplot") ||
                inherits(pl, "ggmatrix") || is.list(pl))
})


## =========================================================
## Section J: CV and Criterion Selection (Book: 02-fect.Rmd)
## =========================================================

test_that("J1: CV selects r minimizing MSPE", {
    skip_on_cran()
    d <- make_panel(N = 50, TT = 20, T0 = 12, Ntr = 10, r = 2, seed = 8110)
    out <- fect(Y ~ D, data = d, index = c("id", "time"),
                method = "ife", r = c(0, 5), se = FALSE, CV = TRUE)
    ## Should have selected an r
    expect_true(!is.null(out$r.cv))
    expect_true(out$r.cv >= 0 && out$r.cv <= 5)
})

test_that("J2: criterion='pc' runs and may select different r", {
    skip_on_cran()
    d <- make_panel(N = 50, TT = 20, T0 = 12, Ntr = 10, r = 2, seed = 8111)
    out_mspe <- fect(Y ~ D, data = d, index = c("id", "time"),
                     method = "ife", r = c(0, 5), se = FALSE,
                     CV = TRUE, criterion = "mspe")
    out_pc <- fect(Y ~ D, data = d, index = c("id", "time"),
                   method = "ife", r = c(0, 5), se = FALSE,
                   CV = TRUE, criterion = "pc")
    ## Both should produce valid r.cv
    expect_true(!is.null(out_mspe$r.cv))
    expect_true(!is.null(out_pc$r.cv))
})


## =========================================================
## Section K: Bundled Dataset Properties (Book: 01-start.Rmd)
## =========================================================

test_that("K1: simdata has 200 units and 30+ time periods", {
    data(simdata, package = "fect")
    expect_equal(length(unique(simdata$id)), 200)
    expect_true(length(unique(simdata$time)) >= 30)
})

test_that("K2: simdata has treatment reversals", {
    data(simdata, package = "fect")
    ## Check that some units switch D from 1 back to 0
    has_reversal <- FALSE
    for (uid in unique(simdata$id)) {
        dvec <- simdata$D[simdata$id == uid]
        diffs <- diff(dvec)
        if (any(diffs == -1)) {
            has_reversal <- TRUE
            break
        }
    }
    expect_true(has_reversal, info = "simdata should contain treatment reversals")
})

test_that("K3: sim_gsynth has 5 treated, 45 control, 30 periods", {
    data(sim_gsynth, package = "fect")
    uid <- unique(sim_gsynth$id)
    expect_equal(length(uid), 50)
    expect_equal(length(unique(sim_gsynth$time)), 30)
    ## 5 treated units (check treatment ever on)
    ever_treated <- tapply(sim_gsynth$D, sim_gsynth$id, max)
    expect_equal(sum(ever_treated == 1), 5)
    expect_equal(sum(ever_treated == 0), 45)
})

test_that("K4: sim_gsynth treatment starts at period 21", {
    data(sim_gsynth, package = "fect")
    ## For treated units, first treatment period should be 21
    ever_treated <- tapply(sim_gsynth$D, sim_gsynth$id, max)
    tr_ids <- names(ever_treated[ever_treated == 1])
    for (uid in tr_ids) {
        subset <- sim_gsynth[sim_gsynth$id == uid, ]
        first_treat <- min(subset$time[subset$D == 1])
        expect_equal(first_treat, 21,
                     info = paste("Unit", uid, "first treatment at", first_treat))
    }
})
