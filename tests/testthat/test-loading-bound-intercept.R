## Bounded factor loadings (loading.bound = "simplex"): the unit intercept.
##
## Behavioural contract under test. With a unit effect in the model
## (force = "unit" or "two-way"), a treated unit's loadings and its intercept
## solve ONE minimisation. So the unit's estimated loadings, weights,
## residuals and treatment effects cannot depend on the unit's outcome level.
## The intercept alone moves, by exactly the constant added to the outcome.
## Without a unit effect (force = "none" or "time") nothing changes.
##
## Groups
##   T1  shift invariance (balanced / unbalanced-T0 / missing-cell panels)
##   T2  noiseless recovery of loading and intercept
##   T3  inactive constraint: bounded fit equals unbounded fit
##   T4  no change for loading.bound = "none" and force "none"/"time"
##   T5  the bound still binds (hull); the intercept is not bounded
##   T6  bootstrap replicates run and carry the invariance
##
## Every panel is built from a fixed seed. Nothing here is tuned to pass.

## ---------------------------------------------------------------------------
## Data generators
## ---------------------------------------------------------------------------

## G1 / G3. Two-factor panel: units 1..Ntr are treated, the rest never are.
## Factor 1 is a standardised linear trend, factor 2 is N(0, 1) noise.
## Loadings and unit effects are N(0, 1) (G1), time effects N(0, 0.5^2),
## errors N(0, 0.3^2), and the treatment effect is +1.
## design = "G3": the treated loadings are the controls' mean loading plus
## N(0, 0.1^2) noise (well inside the controls' hull), unit effects have
## mean 2, and r = 1 uses the trend factor only.
## shape: "balanced"  adoption at period 11 for all treated units;
##        "unbalT0"   adoption at periods 9, 11, 13;
##        "missing"   adoption at 11, outcome NA at (unit 1, periods 3 and 7)
##                    and (unit 2, period 5), rows kept.
lbi_panel <- function(shape = c("balanced", "unbalT0", "missing"),
                      Nco = 20, Ntr = 3, TT = 16, seed = 11,
                      design = c("G1", "G3"), r = 2) {
    shape  <- match.arg(shape)
    design <- match.arg(design)
    N  <- Nco + Ntr
    tt <- seq_len(TT)
    set.seed(seed)
    f1  <- (tt - mean(tt)) / sd(tt)
    f2  <- rnorm(TT)
    Fm  <- cbind(f1, f2)
    Lam <- matrix(rnorm(N * 2), N, 2)
    if (design == "G3") {
        mu <- colMeans(Lam[(Ntr + 1):N, , drop = FALSE])
        Lam[1:Ntr, ] <- matrix(mu, Ntr, 2, byrow = TRUE) +
            matrix(rnorm(Ntr * 2, 0, 0.1), Ntr, 2)
    }
    alpha <- rnorm(N) + if (design == "G3") 2 else 0
    xi    <- rnorm(TT, 0, 0.5)
    err   <- matrix(rnorm(TT * N, 0, 0.3), TT, N)
    nf    <- if (design == "G3" && r == 1) 1 else 2
    T0 <- switch(shape, balanced = c(11, 11, 11), unbalT0 = c(9, 11, 13),
                 missing = c(11, 11, 11))
    D <- matrix(0, TT, N)
    for (i in seq_len(Ntr)) D[tt >= T0[i], i] <- 1
    Y <- Fm[, seq_len(nf), drop = FALSE] %*%
             t(Lam[, seq_len(nf), drop = FALSE]) +
         matrix(alpha, TT, N, byrow = TRUE) + xi + err + D
    if (shape == "missing") {
        Y[c(3, 7), 1] <- NA
        Y[5, 2]       <- NA
    }
    data.frame(id = rep(seq_len(N), each = TT), time = rep(tt, N),
               Y = c(Y), D = c(D))
}

## G2. Noiseless one-factor design. Controls: Y = a_j + b_j F_t exactly.
## Treated unit 9: Y = 1 + 0.5 F_t, treated from period 5, no effect.
## The controls' loadings bracket 0.5 (mean 0.5); `bshift` moves them all.
lbi_noiseless <- function(bshift = 0) {
    Ft <- c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5)
    b  <- c(-0.2, 0.1, 0.3, 0.7, 0.9, 1.2, 0.0, 1.0) + bshift
    a  <- c(-1.5, -0.5, 0.5, 1.5, -1, 1, 0, 0)
    Y  <- cbind(sapply(1:8, function(j) a[j] + b[j] * Ft), 1 + 0.5 * Ft)
    D  <- matrix(0, 6, 9)
    D[5:6, 9] <- 1
    list(data = data.frame(id = rep(1:9, each = 6), time = rep(1:6, 9),
                           Y = c(Y), D = c(D)),
         b = b)
}

## Add constants to chosen units' outcomes; `shifts` is a named list keyed by
## unit id, e.g. list("1" = 10).
lbi_shift <- function(df, shifts) {
    for (u in names(shifts)) {
        k <- df$id == as.integer(u)
        df$Y[k] <- df$Y[k] + shifts[[u]]
    }
    df
}

## Standard gsynth fit. gamma = NULL picks gamma by the package's own
## cross-validation; set.seed(5) first so the folds match across calls.
lbi_fit <- function(df, force, gamma, r = 2, bound = "simplex", ...) {
    set.seed(5)
    suppressMessages(suppressWarnings(fect::fect(
        Y ~ D, data = df, index = c("id", "time"), method = "gsynth",
        force = force, r = r, CV = FALSE, se = FALSE, parallel = FALSE,
        loading.bound = bound, gamma.loading = gamma, ...)))
}

## ---------------------------------------------------------------------------
## Comparison helpers
## ---------------------------------------------------------------------------

## Largest absolute difference; Inf when the NA patterns differ.
lbi_md <- function(a, b) {
    a <- as.numeric(a)
    b <- as.numeric(b)
    if (length(a) != length(b) || any(is.na(a) != is.na(b))) return(Inf)
    if (all(is.na(a))) return(0)
    max(abs(a - b), na.rm = TRUE)
}

expect_small <- function(d, tol, what) {
    testthat::expect(isTRUE(d <= tol),
        sprintf("%s: max abs difference %.3g exceeds tolerance %.0e",
                what, d, tol))
}

## Column of Y.ct that belongs to unit `id` (columns follow colnames(eff)).
lbi_yct <- function(fit, id) {
    fit$Y.ct[, match(as.character(id), colnames(fit$eff))]
}

## Compare a fit on the original data with a fit on data in which unit
## outcomes were shifted by `shifts` (named by unit id). Every quantity that
## should not depend on the level of those outcomes must agree to `tol`;
## the intercept and the imputed outcome move by exactly the constant.
lbi_check_shift <- function(f0, f1, shifts, tol = 1e-8, label = "") {
    p <- function(x) paste(label, x)
    for (s in c("att.avg", "att", "eff", "lambda.tr", "lambda.co",
                "wgt.implied", "loading.proj.resid", "alpha.co", "factor"))
        expect_small(lbi_md(f0[[s]], f1[[s]]), tol, p(s))
    ## gamma is chosen from a grid: it must be exactly the same value.
    expect_equal(f1$gamma.loading, f0$gamma.loading, tolerance = 0,
                 label = p("gamma.loading"))
    ids <- colnames(f0$eff)
    for (u in ids) {
        cst <- if (u %in% names(shifts)) shifts[[u]] else 0
        expect_small(lbi_md(lbi_yct(f1, u) - cst, lbi_yct(f0, u)), tol,
                     p(paste0("Y.ct unit ", u)))
    }
    tr <- rownames(f0$lambda.tr)
    for (k in seq_along(tr)) {
        cst <- if (tr[k] %in% names(shifts)) shifts[[tr[k]]] else 0
        expect_small(abs(f1$alpha.tr[k] - f0$alpha.tr[k] - cst), tol,
                     p(paste0("alpha.tr unit ", tr[k])))
    }
}

## ---------------------------------------------------------------------------
## T1. Shift invariance
## ---------------------------------------------------------------------------

## Property: adding a constant to a treated unit's outcome changes only its
## intercept and its imputed outcome (both by that constant). Loadings,
## weights, residuals, effects, the chosen gamma and every control unit's
## slots stay put. One block per panel shape and gamma rule.

lbi_t1 <- function(shape, gamma, forces = c("unit", "two-way")) {
    d0 <- lbi_panel(shape)
    for (fo in forces) {
        f0 <- lbi_fit(d0, fo, gamma)
        f1 <- lbi_fit(lbi_shift(d0, list("1" = 10)), fo, gamma)
        lbi_check_shift(f0, f1, list("1" = 10), label = paste(shape, fo, "+10:"))
    }
}

test_that("T1: shift invariance, balanced panel, fixed gamma", {
    lbi_t1("balanced", 1)
})
test_that("T1: shift invariance, unbalanced adoption times, fixed gamma", {
    skip_on_cran()
    lbi_t1("unbalT0", 1)
})
test_that("T1: shift invariance, missing cells, fixed gamma", {
    skip_on_cran()
    lbi_t1("missing", 1)
})
test_that("T1: shift invariance, balanced panel, cross-validated gamma", {
    skip_on_cran()
    lbi_t1("balanced", NULL)
})
test_that("T1: shift invariance, unbalanced adoption times, cross-validated gamma", {
    skip_on_cran()
    lbi_t1("unbalT0", NULL)
})
test_that("T1: shift invariance, missing cells, cross-validated gamma", {
    skip_on_cran()
    lbi_t1("missing", NULL)
})

## Property: the same holds for a very large shift, and the intercept is not
## bounded: it rises by the whole shift although the controls' intercepts
## are within a few units of zero.
test_that("T1: shift of 1000 with cross-validated gamma, intercept unbounded", {
    skip_on_cran()
    for (shape in c("balanced", "unbalT0", "missing")) {
        d0 <- lbi_panel(shape)
        for (fo in c("unit", "two-way")) {
            f0 <- lbi_fit(d0, fo, NULL)
            f1 <- lbi_fit(lbi_shift(d0, list("1" = 1000)), fo, NULL)
            lbi_check_shift(f0, f1, list("1" = 1000),
                            label = paste(shape, fo, "+1000:"))
            expect_lt(max(abs(f0$alpha.co)), 10)
            expect_small(abs(f1$alpha.tr[1] - f0$alpha.tr[1] - 1000), 1e-8,
                         paste(shape, fo, "alpha.tr rises by 1000"))
        }
    }
})

## Property: shifting all three treated units at once, each by its own
## constant, leaves every unit's loadings and effects unchanged.
test_that("T1: shifts (10, -3, 40) on the three treated units at once", {
    skip_on_cran()
    d0 <- lbi_panel("balanced")
    sh <- list("1" = 10, "2" = -3, "3" = 40)
    for (fo in c("unit", "two-way")) for (g in list(1, NULL)) {
        f0 <- lbi_fit(d0, fo, g)
        f1 <- lbi_fit(lbi_shift(d0, sh), fo, g)
        lbi_check_shift(f0, f1, sh,
            label = paste("balanced", fo, if (is.null(g)) "CV" else "gamma 1"))
    }
})

## ---------------------------------------------------------------------------
## T2. Noiseless recovery
## ---------------------------------------------------------------------------

## Property: when the data are exactly Y = a + b F, the treated unit's true
## loading lies inside the controls' hull and gamma is huge, the unit's
## loading and intercept are recovered and the post-period prediction error
## is zero. With only 4 pretreatment periods the default min.T0 refuses the
## panel, so min.T0 = 3 is passed.

lbi_noiseless_fit <- function(nl, force) {
    suppressMessages(suppressWarnings(fect::fect(
        Y ~ D, data = nl$data, index = c("id", "time"), method = "gsynth",
        force = force, r = 1, CV = FALSE, se = FALSE, parallel = FALSE,
        loading.bound = "simplex", gamma.loading = 1e6, min.T0 = 3)))
}

test_that("T2: noiseless data, unit effect only: loading, intercept, prediction", {
    nl <- lbi_noiseless()
    ## the true loading 0.5 lies strictly inside the controls' range
    expect_gt(0.5, min(nl$b)); expect_lt(0.5, max(nl$b))
    f <- lbi_noiseless_fit(nl, "unit")
    expect_small(abs(f$att.avg), 1e-4, "att.avg is zero")
    expect_small(lbi_md(lbi_yct(f, 9)[5:6], c(1.75, 2.25)), 1e-4,
                 "post-period prediction")
    expect_small(abs(f$alpha.tr[1] - 1), 1e-4, "intercept is 1")
    ## fitted control loadings are proportional to b (factors are identified
    ## only up to scale); express the treated loading on b's scale
    s <- sum(f$lambda.co * nl$b) / sum(nl$b^2)
    expect_small(lbi_md(f$lambda.co, s * nl$b), 1e-8,
                 "control loadings proportional to b")
    expect_small(abs(f$lambda.tr[1] / s - 0.5), 1e-4, "loading is 0.5")
})

test_that("T2: noiseless data, two-way effects: zero prediction error", {
    for (bs in c(0, 0.1)) {
        nl <- lbi_noiseless(bs)
        f  <- lbi_noiseless_fit(nl, "two-way")
        lab <- paste("controls' loadings shifted by", bs)
        expect_small(abs(f$att.avg), 1e-4, paste(lab, "att.avg is zero"))
        expect_small(lbi_md(lbi_yct(f, 9)[5:6], c(1.75, 2.25)), 1e-4,
                     paste(lab, "post-period prediction"))
        expect_small(abs(f$alpha.tr[1] - 1), 1e-4, paste(lab, "intercept"))
    }
})

## ---------------------------------------------------------------------------
## T3. Inactive constraint
## ---------------------------------------------------------------------------

## True when each row of `Ltr` lies inside the convex hull of the rows of
## `Lco` (r = 1: between smallest and largest; r = 2: inside the polygon).
lbi_in_hull <- function(Lco, Ltr, eps = 1e-10) {
    if (ncol(Lco) == 1) {
        return(all(Ltr >= min(Lco) - eps & Ltr <= max(Lco) + eps))
    }
    h <- grDevices::chull(Lco)
    P <- Lco[h, , drop = FALSE]
    area2 <- sum(P[, 1] * P[c(2:nrow(P), 1), 2] - P[c(2:nrow(P), 1), 1] * P[, 2])
    if (area2 < 0) P <- P[nrow(P):1, , drop = FALSE]   # counter-clockwise
    Q <- P[c(2:nrow(P), 1), , drop = FALSE]
    all(vapply(seq_len(nrow(Ltr)), function(i) {
        q <- Ltr[i, ]
        cr <- (Q[, 1] - P[, 1]) * (q[2] - P[, 2]) -
              (Q[, 2] - P[, 2]) * (q[1] - P[, 1])
        all(cr >= -eps)
    }, logical(1)))
}

## Property: if the treated units' own (unconstrained) loadings lie inside
## the convex hull of the controls' loadings, a very large gamma lets the
## bounded fit reach them, so the bounded fit equals the unbounded fit
## (loadings, intercepts, effects).
test_that("T3: bounded fit equals unbounded fit when the constraint is inactive", {
    skip_on_cran()
    for (r in c(1, 2)) for (shape in c("balanced", "unbalT0", "missing")) {
        d <- lbi_panel(shape, Nco = 25, seed = 21, design = "G3", r = r)
        for (fo in c("unit", "two-way")) {
            lab <- paste("r =", r, shape, fo)
            fu <- lbi_fit(d, fo, NULL, r = r, bound = "none")
            fb <- lbi_fit(d, fo, 1e6, r = r, bound = "simplex")
            ## precondition: the constraint really is inactive
            expect_true(lbi_in_hull(fu$lambda.co, fu$lambda.tr),
                        label = paste(lab, "(precondition: loadings in hull)"))
            expect_small(lbi_md(fb$att.avg, fu$att.avg), 1e-4,
                         paste(lab, "att.avg"))
            expect_small(lbi_md(fb$lambda.tr, fu$lambda.tr), 1e-4,
                         paste(lab, "lambda.tr"))
            expect_small(lbi_md(fb$alpha.tr, fu$alpha.tr), 1e-4,
                         paste(lab, "alpha.tr"))
            expect_small(lbi_md(fb$eff, fu$eff), 1e-3, paste(lab, "eff"))
        }
    }
})

## ---------------------------------------------------------------------------
## T4. No change for the unaffected paths
## ---------------------------------------------------------------------------

## Property: loading.bound = "none", and loading.bound = "simplex" with no
## unit effect (force "none" or "time"), give the results of fect 2.4.7.
## The constants below are att.avg (12 significant digits) from a fect 2.4.7
## build (origin/dev 4282be1), taken 2026-10-16, on lbi_panel(<shape>) with
## seed 11 (G1: Nco = 20, Ntr = 3, TT = 16), r = 2, CV = FALSE, set.seed(5)
## before the fit. A full cross-build comparison of all slots (two R
## processes, tolerance 1e-12) was run once when these values were taken.
test_that("T4: att.avg matches the pinned fect 2.4.7 values", {
    tol <- 1e-8
    pins_simplex <- list(   # loading.bound = "simplex", gamma.loading = 1
        list("balanced", "none", 0.597427322408),
        list("balanced", "time", 0.902076955838),
        list("unbalT0",  "none", 0.638566730003),
        list("unbalT0",  "time", 0.962237652246),
        list("missing",  "none", 0.648439337953),
        list("missing",  "time", 0.923052304804))
    for (p in pins_simplex) {
        f <- lbi_fit(lbi_panel(p[[1]]), p[[2]], 1)
        expect_small(abs(f$att.avg - p[[3]]), tol,
                     paste("simplex", p[[1]], "force", p[[2]]))
    }
    pins_none <- list(      # loading.bound = "none"
        list("balanced", "none", 0.318041057125),
        list("balanced", "unit", 0.414067680502))
    for (p in pins_none) {
        f <- lbi_fit(lbi_panel(p[[1]]), p[[2]], NULL, bound = "none")
        expect_small(abs(f$att.avg - p[[3]]), tol,
                     paste("none", p[[1]], "force", p[[2]]))
    }
})

## ---------------------------------------------------------------------------
## T5. The bound still binds
## ---------------------------------------------------------------------------

## Property, for every force: each treated unit's weights are a point on the
## simplex over the controls, its loadings are those weights applied to the
## controls' loadings, and its factor component lies, period by period,
## between the smallest and largest control factor component. The intercept
## is outside the constraint.
lbi_check_bound <- function(fit, label, Ntr = 3, Nco = 20) {
    W <- fit$wgt.implied
    expect_equal(dim(W), c(Nco, Ntr), label = paste(label, "dim(wgt.implied)"))
    expect_small(lbi_md(colSums(W), rep(1, Ntr)), 1e-6,
                 paste(label, "weights sum to one"))
    expect_gte(min(W), -1e-10)
    for (i in seq_len(Ntr)) {
        expect_small(lbi_md(fit$lambda.tr[i, ],
                            t(fit$lambda.co) %*% W[, i]), 1e-6,
                     paste(label, "lambda.tr = lambda.co' w, unit", i))
    }
    comp_tr <- fit$factor %*% t(fit$lambda.tr)
    comp_co <- fit$factor %*% t(fit$lambda.co)
    lo <- apply(comp_co, 1, min)
    hi <- apply(comp_co, 1, max)
    for (i in seq_len(Ntr)) {
        expect_true(all(comp_tr[, i] >= lo - 1e-6 & comp_tr[, i] <= hi + 1e-6),
                    label = paste(label, "factor component in range, unit", i))
    }
    expect_identical(fit$loading.bound, "simplex")
    expect_true(is.finite(fit$gamma.loading) && fit$gamma.loading > 0,
                label = paste(label, "gamma.loading finite and positive"))
    expect_length(fit$loading.proj.resid, Ntr)
    expect_true(all(is.finite(fit$loading.proj.resid) &
                    fit$loading.proj.resid >= 0),
                label = paste(label, "loading.proj.resid finite, non-negative"))
}

test_that("T5: weights, loadings and factor component stay inside the controls' hull", {
    for (shape in c("balanced", "unbalT0", "missing")) {
        d <- lbi_panel(shape)
        for (fo in c("unit", "two-way", "none", "time")) {
            lbi_check_bound(lbi_fit(d, fo, 1), paste(shape, fo, "gamma 1:"))
        }
    }
})
test_that("T5: same, with cross-validated gamma", {
    skip_on_cran()
    for (shape in c("balanced", "unbalT0", "missing")) {
        d <- lbi_panel(shape)
        for (fo in c("unit", "two-way", "none", "time")) {
            lbi_check_bound(lbi_fit(d, fo, NULL), paste(shape, fo, "CV:"))
        }
    }
})

## Property: with a unit effect, loading.proj.resid[i] is the Euclidean norm
## of the unit's pretreatment residuals from its own (joint) fit, i.e. the
## pretreatment cells of eff (observed, D = 0).
test_that("T5: loading.proj.resid is the norm of the unit's pretreatment residuals", {
    for (shape in c("balanced", "unbalT0", "missing")) {
        d <- lbi_panel(shape)
        for (fo in c("unit", "two-way")) {
            f <- lbi_fit(d, fo, 1)
            for (i in 1:3) {
                dd <- d[d$id == i, ]
                e  <- f$eff[, match(as.character(i), colnames(f$eff))]
                pre <- dd$D == 0 & !is.na(dd$Y)
                expect_small(abs(f$loading.proj.resid[i] - sqrt(sum(e[pre]^2))),
                             1e-6, paste(shape, fo, "unit", i))
            }
        }
    }
})

## ---------------------------------------------------------------------------
## T6. Bootstrap
## ---------------------------------------------------------------------------

## Property: bootstrap and parametric-bootstrap fits with the bound run, drop
## no replicate, give finite standard errors, and the replicate draws are
## (nearly) unchanged by a constant added to a treated unit's outcome. The
## replicate tolerances are set by how each draw is made: 1e-8 for the
## nonparametric draws at a fixed gamma, 1e-4 when gamma is cross-validated
## in each draw, 1e-3 for the parametric draws.
test_that("T6: bootstrap replicates run and are shift-invariant", {
    skip_on_cran()
    d0 <- lbi_panel("balanced")
    d1 <- lbi_shift(d0, list("1" = 10))
    nb <- 15
    run <- function(d, vt, g) {
        set.seed(3)
        suppressMessages(suppressWarnings(fect::fect(
            Y ~ D, data = d, index = c("id", "time"), method = "gsynth",
            force = "two-way", r = 2, CV = FALSE, se = TRUE, nboots = nb,
            seed = 77, parallel = FALSE, vartype = vt,
            loading.bound = "simplex", gamma.loading = g)))
    }
    cfg <- list(list("bootstrap",  1,    1e-8),
                list("bootstrap",  NULL, 1e-4),
                list("parametric", 1,    1e-3),
                list("parametric", NULL, 1e-3))
    for (cf in cfg) {
        lab <- paste(cf[[1]], if (is.null(cf[[2]])) "CV gamma" else "gamma 1")
        f0 <- run(d0, cf[[1]], cf[[2]])
        f1 <- run(d1, cf[[1]], cf[[2]])
        for (f in list(f0, f1)) {
            expect_length(f$att.avg.boot, nb)
            expect_false(anyNA(f$att.avg.boot), label = paste(lab, "no NA draws"))
            expect_true(is.finite(f$est.avg[1, "S.E."]),
                        label = paste(lab, "finite S.E."))
        }
        expect_small(abs(f0$att.avg - f1$att.avg), 1e-8,
                     paste(lab, "point estimate"))
        expect_small(lbi_md(f0$att.avg.boot, f1$att.avg.boot), cf[[3]],
                     paste(lab, "replicate draws"))
        expect_lt(abs(f1$est.avg[1, "S.E."] / f0$est.avg[1, "S.E."] - 1), 1e-3)
    }
})
