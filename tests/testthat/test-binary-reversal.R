## Binary Probit IFE with treatment reversals.

make_binary_reversal_panel <- function(seed = 411) {
    set.seed(seed)
    N <- 24L; TT <- 24L
    d <- expand.grid(time = seq_len(TT), id = seq_len(N))
    on <- as.integer(d$id > 16 & d$time >= 12 + (d$id %% 3))
    ## units 21-24 switch treatment off again after 6 periods
    off <- as.integer(d$id > 20 & d$time >= 18 + (d$id %% 3))
    d$D <- on * (1L - off)
    d$X <- rnorm(nrow(d), sd = 0.3)
    u <- rnorm(N, sd = 0.25); tt <- rnorm(TT, sd = 0.2)
    d$Y <- rbinom(nrow(d), 1, pnorm(d$X + u[d$id] + tt[d$time] + 0.4 * d$D))
    d
}

rev_fit <- function(d, ...) {
    suppressMessages(fect(Y ~ D + X, data = d, index = c("id", "time"),
        binary = TRUE, force = "two-way", QR = TRUE, parallel = FALSE,
        tol = 1e-3, ...))
}

test_that("binary fixed-rank fit runs with reversals and reports switch-off effects", {
    d <- make_binary_reversal_panel()
    D <- matrix(d$D, 24, 24)
    expect_true(any(apply(D, 2, function(v) any(diff(v) < 0))))
    for (r in 0:1) {
        fit <- rev_fit(d, CV = FALSE, r = r, se = FALSE)
        expect_true(is.finite(fit$att.avg))
        expect_true(all(fit$Y.ct.full >= 0 & fit$Y.ct.full <= 1))
        expect_true(length(fit$time.off) >= 1)
        expect_equal(length(fit$att.off), length(fit$time.off))
        ## effects are only defined on observed cells; treated cells enter ATT
        eff <- fit$eff; Dm <- fit$D.dat
        expect_equal(fit$att.avg, mean(eff[Dm == 1 & !is.na(eff)]))
        expect_equal(fit$hasRevs, 1)
    }
})

test_that("binary rolling CV runs with reversals", {
    skip_on_cran()
    d <- make_binary_reversal_panel()
    fit <- suppressWarnings(rev_fit(d, CV = TRUE, r = c(0, 1), se = FALSE,
        cv.method = "rolling", k = 3, cv.nobs = 2, cv.buffer = 1, min.T0 = 5, seed = 3))
    expect_true(fit$r.cv %in% 0:1)
    expect_true(is.matrix(fit$CV.out))
    ## no held-out fold may use a later observed cell of the same unit for training
    TT <- 24L
    for (f in fit$cv.folds) {
        unit <- (f$est.id - 1L) %/% TT + 1L
        t.hold <- (f$est.id - 1L) %% TT + 1L
        II <- fit$I.dat; II[fit$D.dat == 1] <- 0L
        for (u in unique(unit)) {
            last <- max(t.hold[unit == u])
            later <- which(II[, u] == 1L)
            later <- later[later > last]
            if (length(later)) expect_true(all(((u - 1L) * TT + later) %in% f$cv.id))
        }
    }
})

test_that("binary jackknife inference runs with reversals", {
    skip_on_cran()
    d <- make_binary_reversal_panel()
    fit <- rev_fit(d, CV = FALSE, r = 0, se = TRUE, vartype = "jackknife")
    expect_true(is.matrix(fit$est.att))
    expect_true(is.matrix(fit$est.att.off))
    expect_true(all(is.finite(fit$est.avg)))
})
