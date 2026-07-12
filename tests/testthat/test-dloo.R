## Tests for the double leave-one-out (dloo) pre-trend placebo flags on fect().
##
## dloo is a closed-form overlay on the in-sample fit (R/dloo.R): it fills the
## same pre.est.att / pre.att.bound / pre.att.boot slots as `loo`, but without
## re-fitting. These tests verify the overlay reproduces the DEFINITION of the
## estimator (2x2 DiDs of cohort means) to machine precision, that
## `dloo_adjust = TRUE` yields Liu's pre-treatment-average baseline
## ((g-2)/(g-1) rescaling), and that the balanced/staggered guards fire.

skip_on_cran()

## ---- helpers ---------------------------------------------------------------

## balanced, staggered, no-reversal panel
make_panel <- function(TT = 8, cohorts = c(4, 6, Inf), n_per = 25, sd = 0.4,
                       seed = 42) {
    set.seed(seed)
    rows <- list(); id <- 0
    for (g in cohorts) for (u in seq_len(n_per)) {
        id <- id + 1
        y <- rnorm(1) + 0.1 * (1:TT) + rnorm(TT, 0, sd)
        d <- as.integer(!is.infinite(g) & (1:TT) >= g)
        if (!is.infinite(g)) y <- y + d * 2.0
        rows[[length(rows) + 1]] <- data.frame(id = id, time = 1:TT, Y = y, D = d)
    }
    do.call(rbind, rows)
}

## brute-force dloo / average baseline straight from the definition:
## size-weighted aggregate over cohorts of the 2x2-DiD-of-means estimator.
bruteforce_agg <- function(dat, correct = FALSE) {
    firsttreat <- tapply(seq_len(nrow(dat)), dat$id, function(ix) {
        w <- which(dat$D[ix] == 1); if (length(w) == 0) Inf else min(dat$time[ix][w])
    })
    dat$g <- firsttreat[as.character(dat$id)]
    Ybar  <- tapply(dat$Y, list(dat$g, dat$time), mean)
    Nsize <- tapply(dat$id, dat$g, function(x) length(unique(x)))
    gs   <- sort(unique(dat$g[is.finite(dat$g)]))
    allg <- sort(unique(dat$g))
    cell <- list()
    for (g in gs) {
        pre    <- 1:(g - 1)
        later  <- allg[allg > g]
        Nlater <- sum(Nsize[as.character(later)])
        for (t in pre) {
            tp_set <- if (correct) pre else setdiff(pre, t)
            num <- 0
            for (gp in later) for (tp in tp_set) {
                did <- Ybar[as.character(g), as.character(t)] -
                       Ybar[as.character(g), as.character(tp)] -
                       Ybar[as.character(gp), as.character(t)] +
                       Ybar[as.character(gp), as.character(tp)]
                num <- num + Nsize[as.character(gp)] * did
            }
            denom <- (if (correct) (g - 1) else (g - 2)) * Nlater
            cell[[length(cell) + 1]] <- data.frame(
                event.time = t - g + 1, val = num / denom,
                Nsize_g = as.numeric(Nsize[as.character(g)]))
        }
    }
    cell <- do.call(rbind, cell)
    ag <- tapply(seq_len(nrow(cell)), cell$event.time, function(ix)
        sum(cell$val[ix] * cell$Nsize_g[ix]) / sum(cell$Nsize_g[ix]))
    data.frame(event.time = as.numeric(names(ag)), ATT = as.numeric(ag))
}

fit_dloo <- function(dat, correct = FALSE, ...) {
    suppressMessages(fect(
        Y ~ D, data = dat, index = c("id", "time"),
        method = "fe", force = "two-way",
        dloo = TRUE, dloo_adjust = correct,
        parallel = FALSE, nboots = 50, seed = 1, ...))
}

## balanced staggered panel with a subgroup label MIXED within each cohort
make_panel_grp <- function(TT = 8, cohorts = c(4, 6, Inf), n_per = 30,
                           sd = 0.4, seed = 42) {
    set.seed(seed); rows <- list(); id <- 0
    for (g in cohorts) for (u in seq_len(n_per)) {
        id  <- id + 1
        grp <- if (u %% 2 == 0) "A" else "B"
        y   <- rnorm(1) + 0.1 * (1:TT) + (grp == "A") * 0.3 + rnorm(TT, 0, sd)
        d   <- as.integer(!is.infinite(g) & (1:TT) >= g)
        if (!is.infinite(g)) y <- y + d * 2.0
        rows[[length(rows) + 1]] <-
            data.frame(id = id, time = 1:TT, Y = y, D = d, grp = grp)
    }
    do.call(rbind, rows)
}

## Independent brute force of the subgroup-wise dloo: the per-unit raw-Y placebo
## val_j(t) = D_t^j - mean_{t' in B}(D_{t'}^j), D_s^j = Y[s,j] - Ybar_later_g(s),
## averaged within each subgroup. (The FEs cancel in the DiD, so this equals the
## overlay's per-unit contribution; shared later-adopter controls.)
bruteforce_group <- function(dat, correct = FALSE) {
    ft <- tapply(seq_len(nrow(dat)), dat$id, function(ix) {
        w <- which(dat$D[ix] == 1); if (length(w) == 0) Inf else min(dat$time[ix][w])
    })
    dat$g <- ft[as.character(dat$id)]
    Ybar  <- tapply(dat$Y, list(dat$g, dat$time), mean)
    Nsize <- tapply(dat$id, dat$g, function(x) length(unique(x)))
    allg  <- sort(unique(dat$g)); gs <- sort(unique(dat$g[is.finite(dat$g)]))
    ug    <- unique(dat[, c("id", "g", "grp")])
    out <- list()
    for (g in gs) {
        pre <- 1:(g - 1); later <- allg[allg > g]
        if (length(later) == 0) next
        Nl <- sum(Nsize[as.character(later)])
        Ylater <- sapply(pre, function(s)
            sum(Nsize[as.character(later)] *
                Ybar[as.character(later), as.character(s)]) / Nl)
        names(Ylater) <- pre
        for (j in ug$id[ug$g == g]) {
            grp <- ug$grp[ug$id == j]
            Yj  <- dat$Y[dat$id == j & dat$time %in% pre]
            names(Yj) <- dat$time[dat$id == j & dat$time %in% pre]
            Dj  <- Yj[as.character(pre)] - Ylater[as.character(pre)]
            names(Dj) <- pre
            for (t in pre) {
                B <- if (correct) pre else setdiff(pre, t)
                out[[length(out) + 1]] <- data.frame(
                    event.time = t - g + 1, grp = grp,
                    val = as.numeric(Dj[as.character(t)] - mean(Dj[as.character(B)])))
            }
        }
    }
    aggregate(val ~ event.time + grp, do.call(rbind, out), mean)
}

## ---- tests -----------------------------------------------------------------

test_that("dloo overlay reproduces the definition (2x2 DiDs) to machine precision", {
    dat <- make_panel()
    fit <- fit_dloo(dat)
    ref <- bruteforce_agg(dat, correct = FALSE)
    got <- data.frame(event.time = as.numeric(rownames(fit$pre.est.att)),
                      ATT = fit$pre.est.att[, "ATT"])
    m <- merge(got, ref, by = "event.time", suffixes = c(".overlay", ".ref"))
    expect_equal(m$ATT.overlay, m$ATT.ref, tolerance = 1e-10)
})

test_that("dloo_adjust reproduces the pre-treatment-average baseline", {
    dat  <- make_panel()
    fitc <- fit_dloo(dat, correct = TRUE)
    refc <- bruteforce_agg(dat, correct = TRUE)
    got  <- data.frame(event.time = as.numeric(rownames(fitc$pre.est.att)),
                       ATT = fitc$pre.est.att[, "ATT"])
    m <- merge(got, refc, by = "event.time", suffixes = c(".overlay", ".ref"))
    expect_equal(m$ATT.overlay, m$ATT.ref, tolerance = 1e-10)
})

test_that("dloo_adjust == (g-2)/(g-1) rescaling of dloo, cohort by cohort", {
    ## Single-cohort-of-interest panel so the event-time aggregate is one cohort
    ## and the (g-2)/(g-1) factor is unambiguous.
    dat <- make_panel(TT = 8, cohorts = c(5, Inf), n_per = 30)
    fit  <- fit_dloo(dat, correct = FALSE)
    fitc <- fit_dloo(dat, correct = TRUE)
    g <- 5
    a  <- fit$pre.est.att[, "ATT"]
    ac <- fitc$pre.est.att[, "ATT"]
    expect_equal(ac, a * (g - 2) / (g - 1), tolerance = 1e-10)
})

test_that("dloo fills the loo-compatible pre.* slots and inference", {
    dat <- make_panel()
    fit <- fit_dloo(dat)
    expect_true(isTRUE(fit$dloo))
    expect_false(isTRUE(fit$dloo_adjust))
    expect_true(all(c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value",
                      "count.on") %in% colnames(fit$pre.est.att)))
    expect_equal(ncol(fit$pre.att.boot), 50)
    expect_true(all(is.finite(fit$pre.est.att[, "S.E."])))
    expect_false(is.null(fit$dloo.test.out))
    ## rows aligned to the object's own pre-treatment event times (ascending)
    expect_equal(as.numeric(rownames(fit$pre.est.att)),
                 sort(fit$time[fit$time <= 0]))
})

test_that("dloo results plot through the loo view without error", {
    dat <- make_panel()
    fit <- fit_dloo(dat)
    expect_error(plot(fit, loo = TRUE), NA)
})

test_that("dloo inference is the overlay applied inside fect's own bootstrap", {
    ## keep.sims retains the native replicates ONLY for this cross-check; dloo
    ## itself does not need it (the overlay is computed in-loop, no retention).
    dat <- fit_dat <- make_panel()
    fit <- fit_dloo(dat, keep.sims = TRUE)     # nboots = 50
    ## the pre-trend draws are one-per-native-replicate
    expect_false(is.null(fit$eff.boot))
    expect_equal(ncol(fit$pre.att.boot), dim(fit$eff.boot)[3])
    ## each in-loop draw EQUALS the overlay applied post-hoc to that same
    ## case-resampled replicate -- confirms the plumbing feeds the right panel
    b   <- 4
    pseudo <- list(eff = fit$eff.boot[, , b], D.dat = fit$D.boot[, , b],
                   I.dat = fit$I.boot[, , b], rawtime = fit$rawtime,
                   group = fit$group, G = NULL)
    prep <- .dloo_prepare(pseudo, "not-yet-treated")
    ev   <- .dloo_eval(prep, rep(1, prep$N), FALSE)
    pre.term <- sort(fit$time[fit$time <= 0])
    manual   <- ev$att[match(pre.term, ev$event.time)]
    expect_equal(unname(manual), unname(fit$pre.att.boot[, b]), tolerance = 1e-10)
})

test_that("dloo does not force keep.sims (no retained per-replicate panels)", {
    dat <- make_panel()
    fit <- fit_dloo(dat)                        # keep.sims left at default FALSE
    expect_true(is.null(fit$eff.boot))
    expect_true(all(is.finite(fit$pre.est.att[, "S.E."])))
})

test_that("dloo with group reproduces the per-subgroup definition", {
    dat <- make_panel_grp()
    fit <- fit_dloo(dat, group = "grp")
    expect_false(is.null(fit$pre.est.group.output))
    expect_setequal(names(fit$pre.est.group.output), c("A", "B"))
    ref <- bruteforce_group(dat, correct = FALSE)
    for (gname in c("A", "B")) {
        got <- fit$pre.est.group.output[[gname]]$pre.est.att
        gv  <- data.frame(event.time = as.numeric(rownames(got)),
                          ATT = got[, "ATT"])
        m   <- merge(gv, ref[ref$grp == gname, c("event.time", "val")],
                     by = "event.time")
        expect_equal(m$ATT, m$val, tolerance = 1e-10)
    }
})

test_that("dloo_adjust with group reproduces the per-subgroup definition", {
    dat <- make_panel_grp()
    fit <- fit_dloo(dat, correct = TRUE, group = "grp")
    ref <- bruteforce_group(dat, correct = TRUE)
    for (gname in c("A", "B")) {
        got <- fit$pre.est.group.output[[gname]]$pre.est.att
        gv  <- data.frame(event.time = as.numeric(rownames(got)),
                          ATT = got[, "ATT"])
        m   <- merge(gv, ref[ref$grp == gname, c("event.time", "val")],
                     by = "event.time")
        expect_equal(m$ATT, m$val, tolerance = 1e-10)
    }
})

test_that("group leaves the pooled dloo unchanged and partitions it exactly", {
    dat  <- make_panel_grp()
    fitg <- fit_dloo(dat, group = "grp")
    fit0 <- fit_dloo(dat)
    ## pooled point estimate identical with vs without the group breakdown
    expect_equal(fitg$pre.est.att[, "ATT"], fit0$pre.est.att[, "ATT"],
                 tolerance = 1e-12)
    ## count-weighted subgroup ATTs reconstruct the pooled series
    et <- as.numeric(rownames(fitg$pre.est.att))
    recon <- vapply(et, function(e) {
        num <- 0; den <- 0
        for (gname in names(fitg$pre.est.group.output)) {
            go <- fitg$pre.est.group.output[[gname]]$pre.est.att
            r  <- match(e, as.numeric(rownames(go)))
            if (!is.na(r) && !is.na(go[r, "ATT"])) {
                num <- num + go[r, "count.on"] * go[r, "ATT"]
                den <- den + go[r, "count.on"]
            }
        }
        num / den
    }, numeric(1))
    expect_equal(recon, unname(fitg$pre.est.att[, "ATT"]), tolerance = 1e-12)
})

test_that("dloo group slots plot through the loo view per subgroup", {
    dat <- make_panel_grp()
    fit <- fit_dloo(dat, group = "grp")
    expect_error(plot(fit, loo = TRUE, show.group = "A"), NA)
    expect_error(plot(fit, loo = TRUE, show.group = "B"), NA)
})

test_that("dloo requires staggered adoption (no reversal)", {
    dat <- make_panel()
    dat$D[dat$id == 1 & dat$time == 7] <- 0   # induce a reversal
    expect_error(
        suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
                              method = "fe", dloo = TRUE, parallel = FALSE)),
        "staggered")
})

test_that("dloo requires a balanced pre-treatment panel", {
    dat <- make_panel()
    ## drop a pre-treatment cell of a treated (cohort g = 4) unit
    treated4 <- min(dat$id[dat$D == 1 & dat$time == 4])
    dat <- dat[!(dat$id == treated4 & dat$time == 2), ]
    expect_error(
        suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
                              method = "fe", dloo = TRUE, parallel = FALSE,
                              na.rm = TRUE)),
        "balanced")
})

test_that("post-treatment missingness is allowed under dloo", {
    dat <- make_panel()
    ## drop a POST-treatment cell of a cohort g = 4 unit (period 8 > 4)
    treated4 <- min(dat$id[dat$D == 1 & dat$time == 4])
    dat <- dat[!(dat$id == treated4 & dat$time == 8), ]
    expect_error(
        suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
                              method = "fe", dloo = TRUE, parallel = FALSE,
                              na.rm = TRUE)),
        NA)
})

test_that("dloo_adjust requires dloo, and incompatible options error", {
    dat <- make_panel()
    expect_error(
        suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
                              method = "fe", dloo_adjust = TRUE, parallel = FALSE)),
        "dloo_adjust")
    expect_error(
        suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
                              method = "fe", dloo = TRUE, loo = TRUE, parallel = FALSE)),
        "simultaneously")
})

test_that("dloo hard-stops on non-fe methods (ife / mc)", {
    dat <- make_panel()
    expect_error(
        suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
                              method = "ife", dloo = TRUE, parallel = FALSE)),
        "two-way fixed-effects")
    expect_error(
        suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
                              method = "mc", dloo = TRUE, parallel = FALSE)),
        "two-way fixed-effects")
})

test_that("dloo CIs follow ci.method, matching fect's own machinery", {
    dat <- make_panel()
    fn <- fit_dloo(dat)                                   # normal (default)
    fb <- suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
              method = "fe", force = "two-way", dloo = TRUE, parallel = FALSE,
              nboots = 200, seed = 1, ci.method = "basic"))
    att <- fn$pre.est.att[, "ATT"]
    ## normal CIs are symmetric about ATT; basic CIs come from .basic_ci_shifted
    expect_equal((fn$pre.est.att[, "CI.lower"] + fn$pre.est.att[, "CI.upper"]) / 2,
                 att, tolerance = 1e-9)
    cm <- .basic_ci_shifted(fb$pre.est.att[, "ATT"], fb$pre.att.boot, 0.05, FALSE)
    expect_equal(unname(fb$pre.est.att[, "CI.lower"]), unname(cm[, 1]), tolerance = 1e-10)
    expect_equal(unname(fb$pre.est.att[, "CI.upper"]), unname(cm[, 2]), tolerance = 1e-10)
    ## point estimates are identical regardless of CI method
    expect_equal(fb$pre.est.att[, "ATT"], att, tolerance = 1e-12)
})

test_that("plot(fit, dloo = TRUE) renders (distinct dloo view)", {
    dat <- make_panel()
    fit <- fit_dloo(dat)
    expect_error(plot(fit, dloo = TRUE), NA)
})
