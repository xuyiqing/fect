## ---------------------------------------------------------
## #168: method = "cfe" read Q and the gamma groups from the first
## unit, and Z and the kappa groups from the first period, of the
## periods x units arrays that fect() fills from the panel. An absent
## row held 0 there, so one missing row of the first unit, or in the
## first period, changed the model: with Q.type = "linear" the trend
## was 0 in that period for every unit; a unit missing period 1 had
## Z = 0; units missing period 1 shared one kappa group; periods the
## first unit missed shared one gamma group.
##
## The fix builds Z and kappa once per unit, and Q and gamma once per
## period, from the long data before the panel is filled, and stops
## when they vary within a unit or a period.
## ---------------------------------------------------------

drop_rows <- function(data, unit, t) {
    data[!(data$id %in% unit & data$time %in% t), ]
}

cfe_fit <- function(data, ...) {
    fect::fect(Y ~ D, data = data, index = c("id", "time"), method = "cfe",
               force = "two-way", se = FALSE, parallel = FALSE, ...)
}

## simdata: units 101-300, periods 1-35; Z = "L1" with a period-level gamma
cfe_fit_Z <- function(data, ...) {
    fect::fect(Y ~ D + X1 + X2, data = data, index = c("id", "time"),
               method = "cfe", force = "two-way", Z = "L1", gamma = "gamma_t",
               se = FALSE, parallel = FALSE, ...)
}

test_that("Q.type = 'linear': relabeling units 1 and 2 of an unbalanced panel does not change the fit", {

  skip_on_cran()

    data("sim_linear", package = "fect")
    d1 <- drop_rows(sim_linear, 1, 20)          # unit 1 loses pre-treatment period 20
    d.swap <- d1
    d.swap$id <- ifelse(d1$id == 1, 2L, ifelse(d1$id == 2, 1L, d1$id))

    f1 <- cfe_fit(d1, Q.type = "linear")
    fs <- cfe_fit(d.swap, Q.type = "linear")

    expect_equal(fs$att.avg, f1$att.avg, tolerance = 1e-6)
    ## the swapped fit stores unit 2 first
    eff.swap <- fs$eff[, c(2L, 1L, seq_len(ncol(fs$eff))[-(1:2)])]
    expect_equal(unname(eff.swap), unname(f1$eff), tolerance = 1e-6)
})

test_that("Q.type = 'linear': dropping a pre-treatment row of unit 1 moves the ATT as little as dropping it from unit 2", {

  skip_on_cran()

    data("sim_linear", package = "fect")
    f.full <- cfe_fit(sim_linear, Q.type = "linear")
    f.drop1 <- cfe_fit(drop_rows(sim_linear, 1, 20), Q.type = "linear")
    f.drop2 <- cfe_fit(drop_rows(sim_linear, 2, 20), Q.type = "linear")

    ## fad9f72: full 1.006357, drop unit 1 1.095926, drop unit 2 1.006344
    expect_lt(abs(f.drop1$att.avg - f.full$att.avg), 1e-3)
    expect_lt(abs(f.drop2$att.avg - f.full$att.avg), 1e-3)
    expect_lt(abs(f.drop1$att.avg - f.drop2$att.avg), 1e-3)

    ## and it is no longer the fit whose trend is 0 in period 20 for every unit
    d.q <- drop_rows(sim_linear, 1, 20)
    d.q$trend <- d.q$time
    d.q$trend[d.q$time == 20] <- 0
    f.zero <- cfe_fit(d.q, Q = "trend")
    expect_gt(abs(f.drop1$att.avg - f.zero$att.avg), 1e-3)
})

test_that("Z: a unit missing period 1 keeps its Z", {

  skip_on_cran()

    data("simdata", package = "fect")
    s <- simdata
    s$gamma_t <- s$time
    ## unit 101 is the first unit, treated from period 19
    z.t1 <- cfe_fit_Z(drop_rows(s, 101, 1))
    z.t2 <- cfe_fit_Z(drop_rows(s, 101, 2))

    ## fad9f72: full 3.075861, drop period 1 3.132223, drop period 2 3.076764
    expect_lt(abs(z.t1$att.avg - z.t2$att.avg), 5e-3)

    u <- which(z.t1$id == 101)
    post <- which(z.t1$D.dat[, u] == 1)
    expect_lt(max(abs(z.t1$eff[post, u] - z.t2$eff[post, u])), 0.5)

    ## and it is no longer the fit with L1 = 0 for unit 101
    s0 <- drop_rows(s, 101, 1)
    s0$L1[s0$id == 101] <- 0
    z.zero <- cfe_fit_Z(s0)
    expect_gt(abs(z.t1$att.avg - z.zero$att.avg), 1e-3)
})

test_that("kappa groups: two units missing their first-period row keep separate trend loadings", {

  skip_on_cran()

    data("sim_linear", package = "fect")
    dk <- drop_rows(sim_linear, c(1, 2), 1)
    fk <- cfe_fit(dk, Q.type = "linear")

    ## kappa[[1]] is q x N, one column per unit; units in one group share a column
    K <- fk$kappa[[1]]
    expect_equal(dim(K), c(1L, 200L))
    expect_gt(max(abs(K[, 1] - K[, 2])), 1e-8)

    ## an explicit kappa that puts units 1 and 2 in one group is a different model
    dkc <- dk
    dkc$kap <- dkc$id
    dkc$kap[dkc$id == 2] <- 1
    fkc <- cfe_fit(dkc, Q.type = "linear", kappa = "kap")
    expect_equal(max(abs(fkc$kappa[[1]][, 1] - fkc$kappa[[1]][, 2])), 0)
    expect_gt(abs(fk$att.avg - fkc$att.avg), 1e-8)
})

test_that("gamma groups: the first unit missing periods 5 and 10 does not merge those periods' Z coefficients", {

  skip_on_cran()

    data("simdata", package = "fect")
    s <- simdata
    s$gamma_t <- s$time
    zg1 <- cfe_fit_Z(drop_rows(s, 101, c(5, 10)))

    ## gamma[[1]] is T x n_z, one row per period; periods in one group share a row
    G <- zg1$gamma[[1]]
    expect_equal(dim(G), c(35L, 1L))
    expect_gt(max(abs(G[5, ] - G[10, ])), 1e-8)

    ## an explicit gamma that puts periods 5 and 10 in one group is a different model
    g1c <- drop_rows(s, 101, c(5, 10))
    g1c$gamma_t[g1c$time == 10] <- 5
    zg1c <- cfe_fit_Z(g1c)
    expect_equal(max(abs(zg1c$gamma[[1]][5, ] - zg1c$gamma[[1]][10, ])), 0)
    expect_gt(abs(zg1$att.avg - zg1c$att.avg), 1e-8)

    ## the second unit missing the same periods gives a comparable ATT
    ## (fad9f72: 3.081148 for unit 101 vs 3.075868 for unit 102)
    zg2 <- cfe_fit_Z(drop_rows(s, 102, c(5, 10)))
    expect_lt(abs(zg1$att.avg - zg2$att.avg), 2e-3)
})

test_that("Z, Q, gamma and kappa that vary within a unit or a period stop with a message", {

  skip_on_cran()

    data("sim_linear", package = "fect")
    data("simdata", package = "fect")
    s <- simdata
    s$gamma_t <- s$time

    ## Z varies within unit 105
    s.z <- s
    s.z$L1[s.z$id == 105 & s.z$time == 3] <- s.z$L1[s.z$id == 105 & s.z$time == 3] + 1
    expect_error(cfe_fit_Z(s.z), "L1.*unit 105")

    ## Q varies within period 7
    d.q <- sim_linear
    d.q$trend <- d.q$time
    d.q$trend[d.q$id == 3 & d.q$time == 7] <- 0
    expect_error(cfe_fit(d.q, Q = "trend"), "trend.*period 7")

    ## gamma varies within period 6
    s.g <- s
    s.g$gamma_t[s.g$id == 110 & s.g$time == 6] <- 7
    expect_error(cfe_fit_Z(s.g), "gamma_t.*period 6")

    ## kappa varies within unit 4
    d.k <- sim_linear
    d.k$kap <- d.k$id
    d.k$kap[d.k$id == 4 & d.k$time == 2] <- 5
    expect_error(cfe_fit(d.k, Q.type = "linear", kappa = "kap"), "kap.*unit 4")

    ## Z missing in every row of unit 107. With the default na.rm = FALSE
    ## the rows with a missing Z are dropped first and the unit leaves the
    ## panel with fect's usual message, so this needs na.rm = TRUE.
    s.na <- s
    s.na$L1[s.na$id == 107] <- NA
    expect_error(cfe_fit_Z(s.na, na.rm = TRUE), "L1.*unit 107")

    ## a Z that is missing in some rows of a unit is taken from the others
    s.some <- s
    s.some$L1[s.some$id == 107 & s.some$time %in% c(1, 2)] <- NA
    expect_silent(fit.some <- cfe_fit_Z(s.some, na.rm = TRUE))
    expect_true(is.finite(fit.some$att.avg))
})

test_that("fect_mspe() with CFE: a hidden cell in unit 1 scores like a hidden cell in unit 2", {

  skip_on_cran()

    data("sim_linear", package = "fect")
    d.na1 <- sim_linear
    d.na1$Y[d.na1$id == 1 & d.na1$time == 20] <- NA
    d.na2 <- sim_linear
    d.na2$Y[d.na2$id == 2 & d.na2$time == 20] <- NA
    f.na1 <- cfe_fit(d.na1, Q.type = "linear")
    f.na2 <- cfe_fit(d.na2, Q.type = "linear")

    ## the same folds (same seed and dimensions) for both fits; the refits
    ## drop the hidden cells, which on fad9f72 zeroed the trend in period 20
    ## whenever a hidden cell was in unit 1 (MSPE 1.4205 vs 1.3735)
    m1 <- suppressWarnings(fect::fect_mspe(f.na1, seed = 168, cv.method = "rolling", k = 5))
    m2 <- suppressWarnings(fect::fect_mspe(f.na2, seed = 168, cv.method = "rolling", k = 5))
    mspe1 <- m1$scores[["MSPE"]]
    mspe2 <- m2$scores[["MSPE"]]
    expect_lt(abs(mspe1 - mspe2) / mspe2, 0.01)
})

test_that("balanced panels: CFE fits equal fect 2.4.7 at fad9f72", {

  skip_on_cran()

    ## Reference values from the fad9f72 build (before the #168 fix). The
    ## fix changes nothing on a balanced panel, so these hold to rounding.
    check_ref <- function(fit, ref, tol = 1e-10) {
        e <- fit$eff
        expect_equal(fit$att.avg, ref$att.avg, tolerance = tol)
        expect_equal(as.numeric(fit$beta), ref$beta, tolerance = tol)
        expect_equal(head(as.numeric(fit$att), 8), ref$att.head, tolerance = tol)
        expect_equal(tail(as.numeric(fit$att), 8), ref$att.tail, tolerance = tol)
        expect_equal(sum(e, na.rm = TRUE), ref$eff.sum, tolerance = tol)
        expect_equal(sum(e^2, na.rm = TRUE), ref$eff.sumsq, tolerance = tol)
        expect_equal(as.numeric(e[ref$cells]), ref$eff.cells, tolerance = tol)
        expect_equal(fit$mu, ref$mu, tolerance = tol)
        expect_equal(fit$sigma2, ref$sigma2, tolerance = tol)
        if (!is.null(ref$gamma.sum)) {
            expect_equal(dim(fit$gamma[[1]]), ref$gamma.dim)
            expect_equal(sum(sapply(fit$gamma, sum)), ref$gamma.sum, tolerance = tol)
        }
        if (!is.null(ref$kappa.sum)) {
            expect_equal(dim(fit$kappa[[1]]), ref$kappa.dim)
            expect_equal(sum(sapply(fit$kappa, sum)), ref$kappa.sum, tolerance = tol)
        }
    }

    ## sim_linear, Q.type = "linear" (kappa = unit)
    data("sim_linear", package = "fect")
    f.lin <- cfe_fit(sim_linear, Q.type = "linear")
    check_ref(f.lin, list(
        att.avg = 1.0063569035047,
        beta = NA_real_,
        att.head = c(-0.183835693496005, 0.00359135376966502, 0.0560619917779327,
                     -0.0459895537361852, -0.114562704090964, -0.0750284960741648,
                     0.14968868087801, -0.0469812366425354),
        att.tail = c(1.00013199249428, 0.907085902279136, 0.859211735385822,
                     0.839880199683731, 1.02042857318795, 1.00374237565241,
                     1.10102069456438, 1.160975750146),
        eff.sum = 805.084898962493,
        eff.sumsq = 10623.3204349965,
        cells = c(1, 7, 101, 555, 1234, 2500, 3333, 4444, 5001, 10000),
        eff.cells = c(-1.61270007157651, 1.52651628833197, -1.57338015621887,
                      -1.01579280973044, -0.0994957954243603, 1.81973738424735,
                      0.363341199544323, -0.164021136763681, -1.49680153797685,
                      -0.170885926839005),
        kappa.sum = -0.0812435770817471,
        kappa.dim = c(1L, 200L),
        mu = 0.642453271170292,
        sigma2 = 0.997109733055106
    ))

    ## sim_trend, Q.type = "bspline" (kappa = unit). This fit reaches
    ## max.iteration = 5000 on fad9f72 as well; the warnings are that.
    data("sim_trend", package = "fect")
    f.bs <- suppressWarnings(cfe_fit(sim_trend, Q.type = "bspline"))
    check_ref(f.bs, list(
        att.avg = 0.310195845422441,
        beta = NA_real_,
        att.head = c(-0.0929360402146247, 0.0730352991726366, 0.105984560584928,
                     -0.0136984865535476, -0.0980568705166484, -0.0725041077952142,
                     0.139994333172511, -0.0671710214019835),
        att.tail = c(0.652666024692629, 0.458918055005335, 0.295323213530783,
                     0.144199727447466, 0.175834728581284, -0.00793453074883486,
                     -0.0969558148985681, -0.243562808645487),
        eff.sum = 248.16689085489,
        eff.sumsq = 158838.255048992,
        cells = c(1, 7, 101, 555, 1234, 2500, 3333, 4444, 5001, 10000),
        eff.cells = c(-0.102805522178803, 0.852037973092854, -0.976288692386803,
                      -0.770714328547949, 0.143793860012085, -21.1861265950661,
                      0.520571360170945, -0.175102572540672, -0.378049588006472,
                      -0.550182490706395),
        kappa.sum = -5.04636707638765,
        kappa.dim = c(5L, 200L),
        mu = 0.305693592397544,
        sigma2 = 0.902347664251554
    ))

    ## simdata, Z = "L1" (gamma = period)
    data("simdata", package = "fect")
    f.z <- fect::fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                      method = "cfe", force = "two-way", Z = "L1",
                      se = FALSE, parallel = FALSE)
    check_ref(f.z, list(
        att.avg = 3.07586120630162,
        beta = c(1.0246398775412, 2.9436414420185),
        att.head = c(-2.07355905213599, -0.284089376432536, -0.729116858884357,
                     -0.884553780481902, -1.12167473648319, -1.41526881360137,
                     -1.14630260083652, -0.996191864359285),
        att.tail = c(12.4121712865802, 9.63287716125647, 12.9780776522333,
                     8.8958394253501, 10.9707732761286, 12.6592651579873,
                     9.37764854237341, 9.10230042374602),
        eff.sum = 4623.08033271483,
        eff.sumsq = 81728.6358836961,
        cells = c(1, 7, 101, 555, 1234, 2500, 3333, 4444, 5001, 7000),
        eff.cells = c(-3.89234750725061, 2.51367184458521, 9.59112414478024,
                      3.28998404092541, -0.700198622198422, 3.04987865093685,
                      0.765829543561081, 6.23747648742699, -1.47105734151522,
                      1.86142496538101),
        gamma.sum = 16.3761223493663,
        gamma.dim = c(35L, 1L),
        mu = 13.5816094749365,
        sigma2 = 7.88527351032461
    ))
})
