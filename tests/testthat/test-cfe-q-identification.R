## ---------------------------------------------------------
## #173: with method = "cfe", a Q column whose mass lies almost
## entirely in a unit's treated (imputed) periods leaves that unit's
## loading (kappa) identified only from imputed cells. The EM then
## contracts along that direction at the missing-information rate,
## 1 - (share of the column observed), so on sim_trend with the
## default B-spline basis (knots at the tertiles of the long time
## vector, 17 and 34; the fifth basis function is supported on periods
## 35-50 and treatment starts in 41) the rate is 0.9996: tol = 1e-5
## needs about 28,000 iterations and the answer is an extrapolation.
##
## The fix is a diagnostic before the first iteration: for each kappa
## group, the largest eigenvalue of I - (Q'Q)^{-1} Q_obs' Q_obs, where
## Q_obs' Q_obs is averaged over the group's units on their observed
## untreated periods. Above 0.9 fect warns once, naming the column,
## the rate and the iteration count, and the max.iteration warning
## repeats the cause. No number changes.
## ---------------------------------------------------------

cfe_fit <- function(data, ...) {
    fect::fect(Y ~ D, data = data, index = c("id", "time"), method = "cfe",
               force = "two-way", se = FALSE, parallel = FALSE, ...)
}

## every warning a call raises, as a character vector
collect_warnings <- function(expr) {
    msgs <- character(0)
    val <- withCallingHandlers(expr, warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    list(value = val, warnings = msgs)
}

## sim_trend as the arrays fect_cfe() receives: 50 periods x 200 units,
## Q the default B-spline basis (df = 5, degree = 3, built on the long
## time vector as fect() does, so the knots are 17 and 34), one kappa
## group per unit, II = 1 on untreated cells
sim_trend_arrays <- function() {
    data("sim_trend", package = "fect")
    d <- sim_trend[order(sim_trend$id, sim_trend$time), ]
    TT <- length(unique(d$time))
    N <- length(unique(d$id))
    Q <- unique(splines::bs(d$time, df = 5, degree = 3, intercept = FALSE))
    X.Q <- array(0, c(TT, N, ncol(Q)),
                 dimnames = list(NULL, NULL, paste0("time.bs", 1:5)))
    for (j in seq_len(ncol(Q))) X.Q[, , j] <- Q[, j]
    X.kappa <- array(rep(seq_len(N), each = TT), c(TT, N, 1),
                     dimnames = list(NULL, sort(unique(d$id)), "id"))
    II <- 1 - matrix(d$D, TT, N)
    list(X.Q = X.Q, X.kappa = X.kappa, kappaQ.id = list(1:5), II = II)
}

test_that(".cfe_q_identification(): a fully observed unit has rate 0, an unobserved column rate 1", {

  skip_on_cran()

    q.id <- fect:::.cfe_q_identification
    ## one column, five periods; unit 1 observed everywhere, unit 2 misses
    ## periods 4-5, where the column has all its mass
    Q <- c(0, 0, 0, 1, 1)
    X.Q <- array(Q, c(5, 2, 1), dimnames = list(NULL, NULL, "step"))
    X.kappa <- array(rep(1:2, each = 5), c(5, 2, 1),
                     dimnames = list(NULL, c("a", "b"), "id"))
    II.full <- matrix(1, 5, 2)
    d0 <- q.id(X.Q, X.kappa, list(1L), II.full, force = 0)
    expect_equal(d0$rate, 0, tolerance = 1e-12)

    II.miss <- II.full
    II.miss[4:5, 2] <- 0
    d1 <- q.id(X.Q, X.kappa, list(1L), II.miss, force = 0)
    expect_equal(d1$rate, 1, tolerance = 1e-12)
    expect_equal(d1$units, "b")
    expect_equal(unname(d1$observed["step"]), 0)
    expect_equal(d1$flagged, "step")

    ## pooling the two units in one kappa group halves the missing mass
    X.kappa1 <- array(1, c(5, 2, 1), dimnames = list(NULL, c("a", "b"), "id"))
    d2 <- q.id(X.Q, X.kappa1, list(1L), II.miss, force = 0)
    expect_equal(d2$rate, 0.5, tolerance = 1e-12)
    expect_equal(d2$units, c("a", "b"))

    ## no Q: nothing to report
    expect_null(q.id(array(0, c(5, 2, 0)), array(0, c(5, 2, 0)), list(),
                     II.full, force = 3))
})

test_that(".cfe_q_identification(): sim_trend with the default B-spline has rate 0.9996 on time.bs5", {

  skip_on_cran()

    a <- sim_trend_arrays()
    d <- fect:::.cfe_q_identification(a$X.Q, a$X.kappa, a$kappaQ.id, a$II,
                                      force = 3)
    ## 0.999595 with the unit intercept (0.99952 with knots at 17.33 and
    ## 33.67, the basis of the diagnosis in #173, built on the 50 unique
    ## times)
    expect_equal(d$rate, 0.999595, tolerance = 1e-5)
    expect_equal(d$flagged, "time.bs5")
    expect_lt(d$observed["time.bs5"], 0.01)
    expect_gt(d$observed["time.bs1"], 0.99)
    ## the worst group is a treated unit (1-80 are treated in 41-50)
    expect_true(all(as.integer(d$units) %in% 1:80))

    msg <- fect:::.cfe_q_identification_message(d, tol = 1e-5)
    expect_match(msg, "time.bs5")
    expect_match(msg, "0.9996")
    expect_match(msg, "28,000 iterations")
    expect_match(msg, "Q.type = c\\(\"linear\", \"quadratic\"\\)")
})

test_that("sim_trend, Q.type = 'bspline': fect warns before iterating and the max.iteration warning names the cause", {

  skip_on_cran()

    data("sim_trend", package = "fect")
    r <- collect_warnings(cfe_fit(sim_trend, Q.type = "bspline",
                                  max.iteration = 50))
    w <- r$warnings
    pre <- grepl("nearly unidentified", w) & grepl("time.bs5", w) &
        grepl("0.9996", w) & grepl("28,000 iterations", w)
    expect_equal(sum(pre), 1L)
    ## the diagnostic comes before the convergence warnings
    expect_equal(which(pre), 1L)
    mi <- grepl("CFE optimization did not converge within 50", w)
    expect_equal(sum(mi), 1L)
    expect_match(w[mi], "time.bs5")
    expect_match(w[mi], "identified from the pre-treatment periods")
    ## the other two are the existing FE-only and EM convergence notes
    expect_equal(length(w), 4L)
    expect_true(any(grepl("FE-only", w)))
    expect_true(any(grepl("EM did not converge", w)))

    ## the numbers are those of the base install (0e100bd), 51 iterations
    f <- r$value
    expect_equal(f$niter, 51)
    expect_equal(f$att.avg, 0.866173599963474, tolerance = 1e-10)
    expect_equal(sum(f$eff, na.rm = TRUE), 691.980238847933, tolerance = 1e-10)
    expect_equal(sum(unlist(f$kappa)), -30.4409374647323, tolerance = 1e-10)
})

test_that("sim_trend, Q.type = 'linear' and c('linear', 'quadratic'): no identification warning, converged numbers unchanged", {

  skip_on_cran()

    data("sim_trend", package = "fect")
    r <- collect_warnings(cfe_fit(sim_trend, Q.type = "linear"))
    expect_length(r$warnings, 0L)
    f <- r$value
    expect_equal(f$niter, 77)
    expect_equal(f$att.avg, 0.605976465879664, tolerance = 1e-10)
    expect_equal(sum(f$eff, na.rm = TRUE), 484.781125219571, tolerance = 1e-10)
    expect_equal(f$sigma2, 1.00087621503388, tolerance = 1e-10)

    ## the book's replacement example (05-cfe.Rmd)
    r2 <- collect_warnings(cfe_fit(sim_trend, Q.type = c("linear", "quadratic")))
    expect_length(r2$warnings, 0L)
    f2 <- r2$value
    expect_equal(f2$niter, 107)
    expect_equal(f2$att.avg, 1.09503820147888, tolerance = 1e-10)
    expect_equal(sum(f2$eff, na.rm = TRUE), 876.026980352639, tolerance = 1e-10)
    expect_equal(sum(f2$eff^2, na.rm = TRUE), 10800.1660232387, tolerance = 1e-10)
    expect_equal(sum(unlist(f2$kappa)), 0.51740180206012, tolerance = 1e-10)
    expect_equal(f2$mu, 0.257350756700497, tolerance = 1e-10)
})

test_that("balanced panels that converged before are unchanged to 1e-10 (sim_linear linear; simdata with Z)", {

  skip_on_cran()

    data("sim_linear", package = "fect")
    r <- collect_warnings(cfe_fit(sim_linear, Q.type = "linear"))
    expect_length(r$warnings, 0L)
    f <- r$value
    expect_equal(f$niter, 75)
    expect_equal(f$att.avg, 1.0063569035047, tolerance = 1e-10)
    expect_equal(sum(f$eff, na.rm = TRUE), 805.084898962493, tolerance = 1e-10)
    expect_equal(sum(f$eff^2, na.rm = TRUE), 10623.3204349965, tolerance = 1e-10)
    expect_equal(sum(unlist(f$kappa)), -0.0812435770817471, tolerance = 1e-10)
    expect_equal(f$mu, 0.642453271170292, tolerance = 1e-10)
    expect_equal(f$sigma2, 0.997109733055106, tolerance = 1e-10)

    data("simdata", package = "fect")
    r2 <- collect_warnings(
        fect::fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                   method = "cfe", force = "two-way", Z = "L1",
                   se = FALSE, parallel = FALSE))
    expect_length(r2$warnings, 0L)
    f2 <- r2$value
    expect_equal(f2$niter, 38)
    expect_equal(f2$att.avg, 3.07586120630162, tolerance = 1e-10)
    expect_equal(as.vector(f2$beta), c(1.0246398775412, 2.9436414420185),
                 tolerance = 1e-10)
    expect_equal(sum(f2$eff, na.rm = TRUE), 4623.08033271483, tolerance = 1e-10)
    expect_equal(sum(f2$eff^2, na.rm = TRUE), 81728.6358836961, tolerance = 1e-10)
    expect_equal(sum(unlist(f2$gamma)), 16.3761223493663, tolerance = 1e-10)
    expect_equal(f2$mu, 13.5816094749365, tolerance = 1e-10)
    expect_equal(f2$sigma2, 7.88527351032461, tolerance = 1e-10)
})
