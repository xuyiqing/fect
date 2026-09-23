## Bai-Ng style PC criterion for the Probit IFE model (fect 1.1.x "IC").

make_pc_panel <- function(seed = 410) {
    set.seed(seed)
    N <- 20L; TT <- 22L
    d <- expand.grid(time = seq_len(TT), id = seq_len(N))
    d$D <- as.integer(d$id > 15 & d$time >= 15 + (d$id %% 3))
    d$X <- rnorm(nrow(d), sd = 0.3)
    u <- rnorm(N, sd = 0.25); t <- rnorm(TT, sd = 0.2)
    d$Y <- rbinom(nrow(d), 1, pnorm(d$X + u[d$id] + t[d$time] + 0.3 * d$D))
    d
}
pc_fit <- function(d, ...) suppressMessages(fect(Y ~ D + X, data = d,
    index = c("id", "time"), binary = TRUE, force = "two-way", QR = TRUE,
    parallel = FALSE, tol = 1e-3, ...))

test_that("fixed-rank binary fit returns PC matching the Bai-Ng formula", {
    d <- make_pc_panel()
    N <- 20; TT <- 22
    for (r in 0:1) {
        fit <- pc_fit(d, CV = FALSE, r = r, se = FALSE)
        expect_true(is.finite(fit$PC))
        expect_equal(fit$PC,
            r * (N + TT) / (N * TT) * log(N * TT / (N + TT)) - 2 * fit$loglikelihood)
    }
})

test_that("binary CV.out carries PC and it does not depend on the split", {
    skip_on_cran()
    d <- make_pc_panel()
    a <- suppressWarnings(pc_fit(d, CV = TRUE, r = c(0, 1), se = FALSE, k = 3,
        cv.method = "rolling", cv.nobs = 2, cv.buffer = 1, min.T0 = 5, seed = 1))
    b <- suppressWarnings(pc_fit(d, CV = TRUE, r = c(0, 1), se = FALSE, k = 2,
        cv.method = "rolling", cv.nobs = 3, cv.buffer = 1, min.T0 = 5, seed = 99))
    expect_true(all(c("r", "IC", "PC", "Log-likelihood", "MSPE") %in% colnames(a$CV.out)))
    expect_equal(a$CV.out[, "PC"], b$CV.out[, "PC"])
    expect_equal(a$CV.out[, "Log-likelihood"], b$CV.out[, "Log-likelihood"])
    ## consistent with the fixed-rank fit at each rank
    for (i in seq_len(nrow(a$CV.out))) {
        fixed <- pc_fit(d, CV = FALSE, r = a$CV.out[i, "r"], se = FALSE)
        expect_equal(unname(a$CV.out[i, "PC"]), fixed$PC, tolerance = 1e-6)
    }
    ## selected-rank fit exposes the same PC
    expect_equal(a$PC, unname(a$CV.out[a$CV.out[, "r"] == a$r.cv, "PC"]), tolerance = 1e-6)
})
