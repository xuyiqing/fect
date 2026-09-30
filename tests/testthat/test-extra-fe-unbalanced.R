## ---------------------------------------------------------
## #171: the extra fixed-effect labels (index[3:], or group.fe) were
## read from the periods x units array that fect() fills from the
## panel, where an absent row held 0. In the C++ every cell belongs to
## the group its label names, so the absent cells of an unbalanced
## panel formed a phantom group with label 0; the never-treated path
## read each unit's label from period 1, so a treated unit missing
## period 1 stopped with "levels in treated units not found in
## controls: 0".
##
## The fix reads the labels per cell from the long data, with NA at
## absent cells; the C++ takes each group's mean over observed cells
## only; the never-treated path reads each unit's labels from its
## observed cells.
## ---------------------------------------------------------

drop_rows <- function(data, unit, t) {
    data[!(data$id %in% unit & data$time %in% t), ]
}

## sim_linear: 200 units, 50 periods, 80 treated units; grp = id %% 5
sim_grp <- function() {
    data("sim_linear", package = "fect", envir = environment())
    d <- sim_linear
    d$grp <- d$id %% 5
    d
}

grp_fit <- function(data, tcf, force = "two-way", ...) {
    fect::fect(Y ~ D, data = data, index = c("id", "time", "grp"),
               method = "cfe", force = force, time.component.from = tcf,
               se = FALSE, parallel = FALSE, ...)
}

## the mean post-treatment effect of one unit
unit_post_mean <- function(fit, unit) {
    k <- which(fit$id == unit)
    post <- which(fit$D.dat[, k] == 1)
    mean(fit$eff[post, k])
}

test_that("never-treated path: a treated unit missing period 1 runs and moves the ATT as little as one missing period 2", {

  skip_on_cran()

    d <- sim_grp()
    u <- sort(unique(d$id[d$D == 1]))[3]

    f.full <- suppressMessages(grp_fit(d, "nevertreated"))
    ## 0e100bd: "Extra fixed effect dimension 1 has levels in treated
    ## units not found in controls: 0"
    f.t1 <- suppressMessages(grp_fit(drop_rows(d, u, 1), "nevertreated"))
    f.t2 <- suppressMessages(grp_fit(drop_rows(d, u, 2), "nevertreated"))

    ## 0e100bd: full 2.289935, drop period 2 2.289802
    expect_lt(abs(f.t1$att.avg - f.t2$att.avg), 1e-3)
    expect_lt(abs(f.t1$att.avg - f.full$att.avg), 1e-3)
    expect_false(anyNA(f.t1$eff[f.t1$I.dat == 1]))
})

test_that("not-yet-treated path: a grouping nested in the units is absorbed by the unit effects on an unbalanced panel", {

  skip_on_cran()

    ## (the issue's 0.054 shift of unit 3's mean post effect after a
    ## period-1 drop is the same without any extra fixed effect: the
    ## unit's period-1 residual under the two-way model is -2.09 over 40
    ## pre-treatment periods)
    d <- sim_grp()
    u <- sort(unique(d$id[d$D == 1]))[3]
    d1 <- drop_rows(d, u, 1)

    f.grp <- suppressMessages(grp_fit(d1, "notyettreated"))
    f.none <- suppressMessages(fect::fect(
        Y ~ D, data = d1, index = c("id", "time"), method = "cfe",
        force = "two-way", time.component.from = "notyettreated",
        se = FALSE, parallel = FALSE))
    expect_equal(f.grp$att.avg, f.none$att.avg, tolerance = 1e-8)
    expect_equal(unname(f.grp$eff), unname(f.none$eff), tolerance = 1e-8)
    expect_false(anyNA(f.grp$eff[f.grp$I.dat == 1]))
})

test_that("not-yet-treated path: with time effects only, a grouping nested in the units gives the least squares fit on an unbalanced panel", {

  skip_on_cran()

    ## a group effect the unit effects do not absorb; the fit on the
    ## untreated cells must be the three-way least squares fit, and the
    ## counterfactual its prediction
    d <- sim_grp()
    d$Y <- d$Y + c(0, 2, -1, 3, 1.5)[d$grp + 1]
    tr <- sort(unique(d$id[d$D == 1]))
    d1 <- drop_rows(d, tr[3], 1)
    d20 <- drop_rows(d, tr[1:20], 1:3)

    ls_att <- function(data) {
        ctrl <- data[data$D == 0, ]
        m <- fixest::feols(Y ~ 1 | time + grp, data = ctrl)
        trd <- data[data$D == 1, ]
        mean(trd$Y - predict(m, newdata = trd))
    }
    for (data in list(d, d1, d20)) {
        f <- suppressMessages(grp_fit(data, "notyettreated", force = "time"))
        ## 0e100bd: d1 3.42265 against 3.42298 (the absent cell, in a
        ## phantom group, was imputed without its group effect)
        expect_equal(f$att.avg, ls_att(data), tolerance = 1e-6)
    }
})

test_that("relabeling units of the unbalanced panel does not change the fit (both paths)", {

  skip_on_cran()

    d <- sim_grp()
    tr <- sort(unique(d$id[d$D == 1]))
    u <- tr[3]
    d1 <- drop_rows(d, u, 1)
    ## give the unit missing period 1 the smallest id, and unit 1 its id
    d.swap <- d1
    d.swap$id <- ifelse(d1$id == u, 1L, ifelse(d1$id == 1L, u, d1$id))

    for (tcf in c("nevertreated", "notyettreated")) {
        f1 <- suppressMessages(grp_fit(d1, tcf))
        fs <- suppressMessages(grp_fit(d.swap, tcf))
        expect_equal(fs$att.avg, f1$att.avg, tolerance = 1e-8)
        expect_equal(unname(fs$att), unname(f1$att), tolerance = 1e-8)
        ## the swapped fit stores the unit under the other id
        k1 <- which(f1$id == u)
        ks <- which(fs$id == 1L)
        expect_equal(unname(fs$eff[, ks]), unname(f1$eff[, k1]),
                     tolerance = 1e-8)
    }
})

test_that("a label that varies within a unit is read per cell", {

  skip_on_cran()

    ## a cell-level grouping: region x half of the panel
    d <- sim_grp()
    d$cell <- (d$id %% 3) * 10 + as.integer(d$time > 25)
    u <- sort(unique(d$id[d$D == 1]))[3]
    d1 <- drop_rows(d, u, 1)

    f <- suppressMessages(fect::fect(
        Y ~ D, data = d1, index = c("id", "time", "cell"), method = "cfe",
        force = "two-way", time.component.from = "notyettreated",
        se = FALSE, parallel = FALSE))
    expect_false(anyNA(f$eff[f$I.dat == 1]))

    ## the same labels under other names give the same fit
    d2 <- d1
    d2$cell <- match(d2$cell, sort(unique(d2$cell))) * 100
    f2 <- suppressMessages(fect::fect(
        Y ~ D, data = d2, index = c("id", "time", "cell"), method = "cfe",
        force = "two-way", time.component.from = "notyettreated",
        se = FALSE, parallel = FALSE))
    expect_equal(f2$att.avg, f$att.avg, tolerance = 1e-10)
    expect_equal(unname(f2$eff), unname(f$eff), tolerance = 1e-10)
})

test_that("never-treated path: a label that varies within a unit is applied per cell", {

  skip_on_cran()

    ## a cell-level grouping with an effect that changes at period 26;
    ## the counterfactual of a treated cell must carry the effect of
    ## the cell's own label, not of the unit's period-1 label
    d <- sim_grp()
    d$cell <- (d$id %% 3) * 10 + as.integer(d$time > 25)
    ceff <- c("0" = 0, "1" = 4, "10" = -3, "11" = 2, "20" = 1, "21" = -2)
    d$Y <- d$Y + ceff[as.character(d$cell)]
    tr <- sort(unique(d$id[d$D == 1]))

    for (data in list(d, drop_rows(d, tr[3], 1))) {
        f <- suppressMessages(fect::fect(
            Y ~ D, data = data, index = c("id", "time", "cell"),
            method = "cfe", force = "time",
            time.component.from = "nevertreated", se = FALSE,
            parallel = FALSE))
        ctrl <- data[data$D == 0, ]
        trd <- data[data$D == 1, ]
        m <- fixest::feols(Y ~ 1 | time + cell, data = ctrl)
        eff.ls <- trd$Y - predict(m, newdata = trd)
        eff.f <- f$eff[cbind(match(trd$time, as.numeric(rownames(f$eff))),
                             match(trd$id, f$id))]
        ## 0e100bd, balanced: ATT 3.4437 against 3.4565, and cells off by
        ## up to 2.5 (the period-1 label was applied to every period);
        ## unbalanced: the error "levels in treated units not found in
        ## controls: 0"
        expect_lt(abs(f$att.avg - mean(eff.ls)), 1e-3)
        expect_lt(max(abs(eff.f - eff.ls)), 0.1)
    }
})

test_that("balanced panels with extra fixed effects are unchanged (0e100bd literals)", {

  skip_on_cran()

    d <- sim_grp()
    f.nt <- suppressMessages(grp_fit(d, "nevertreated"))
    f.nyt <- suppressMessages(grp_fit(d, "notyettreated"))
    expect_equal(f.nt$att.avg, 2.28993476888885, tolerance = 1e-10)
    expect_equal(unname(f.nt$att[1:3]),
                 c(-1.30758257364332, -0.943861080709907, -0.80506693653883),
                 tolerance = 1e-10)
    expect_equal(f.nyt$att.avg, 2.28993476888885, tolerance = 1e-10)
    expect_equal(unname(f.nyt$att[1:3]),
                 c(-0.784549544185994, -0.566316648425944, -0.483040161923299),
                 tolerance = 1e-10)

    ## simdata with covariates and one factor: two extra dimensions, one
    ## varying within a unit (not-yet-treated)
    data("simdata", package = "fect", envir = environment())
    s <- simdata
    s$grp <- s$id %% 4
    s$grp2 <- (s$id %% 3) * 10 + (s$time > 20)
    f.s1 <- suppressMessages(fect::fect(
        Y ~ D + X1 + X2, data = s, index = c("id", "time", "grp", "grp2"),
        method = "cfe", force = "two-way", r = 1,
        time.component.from = "notyettreated", se = FALSE, parallel = FALSE))
    expect_equal(f.s1$att.avg, 3.36344785399659, tolerance = 1e-10)
    expect_equal(unname(c(f.s1$beta)), c(1.04201840350418, 2.97667934970522),
                 tolerance = 1e-10)

    ## sim_linear with one factor: two unit-level dimensions (never-treated)
    d$grp3 <- d$id %% 3
    f.nt2 <- suppressMessages(fect::fect(
        Y ~ D, data = d, index = c("id", "time", "grp", "grp3"),
        method = "cfe", force = "two-way", r = 1,
        time.component.from = "nevertreated", se = FALSE, parallel = FALSE))
    expect_equal(f.nt2$att.avg, 2.1700383622831, tolerance = 1e-10)
    expect_equal(unname(f.nt2$att[1:3]),
                 c(-1.32503181405856, -0.842810931568849, -0.607219244900253),
                 tolerance = 1e-10)
})
