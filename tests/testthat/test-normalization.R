## ---------------------------------------------------------
## Normalization regression test
## Verifies that sigma2 (est.fect) is consistent with and
## without normalize=TRUE (the sigma2 normalization bug fix)
## ---------------------------------------------------------

test_that("sigma2 is consistent with and without normalization (FE)", {

  skip_on_cran()

    set.seed(6001)

    N <- 30
    TT <- 15
    T0 <- 10
    Ntr <- 10
    tau <- 3.0

    alpha_i <- rnorm(N, 0, 1)
    xi_t <- rnorm(TT, 0, 0.5)

    Y_vec <- numeric(N * TT)
    D_vec <- integer(N * TT)
    id_vec <- integer(N * TT)
    time_vec <- integer(N * TT)

    idx <- 1
    for (i in 1:N) {
        for (t in 1:TT) {
            treated <- (i <= Ntr) && (t > T0)
            D_vec[idx] <- as.integer(treated)
            eps <- rnorm(1, 0, 1)
            Y_vec[idx] <- alpha_i[i] + xi_t[t] + tau * D_vec[idx] + eps
            id_vec[idx] <- i
            time_vec[idx] <- t
            idx <- idx + 1
        }
    }

    simdf <- data.frame(
        id = id_vec,
        time = time_vec,
        Y = Y_vec,
        D = D_vec
    )

    out_raw <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "fe",
        force = "two-way",
        normalize = FALSE,
        se = FALSE,
        parallel = FALSE
    ))

    out_norm <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "fe",
        force = "two-way",
        normalize = TRUE,
        se = FALSE,
        parallel = FALSE
    ))

    ## ATT should be essentially the same
    expect_equal(out_raw$att.avg, out_norm$att.avg, tolerance = 0.01,
        label = "ATT should match with/without normalization")

    ## sigma2.fect should be close (within 10% relative error)
    ## This is the regression test for the normalization bug fix
    if (!is.null(out_raw$sigma2.fect) && !is.null(out_norm$sigma2.fect)) {
        rel_diff <- abs(out_raw$sigma2.fect - out_norm$sigma2.fect) /
            max(abs(out_raw$sigma2.fect), 1e-10)
        expect_lt(rel_diff, 0.1,
            label = paste0("sigma2.fect relative difference = ",
                           round(rel_diff, 6),
                           " (raw=", round(out_raw$sigma2.fect, 4),
                           ", norm=", round(out_norm$sigma2.fect, 4), ")"))
    }
})

test_that("sigma2 is consistent with and without normalization (IFE)", {

  skip_on_cran()

    set.seed(6002)

    N <- 30
    TT <- 15
    T0 <- 10
    Ntr <- 10
    tau <- 2.0

    alpha_i <- rnorm(N, 0, 1)
    xi_t <- rnorm(TT, 0, 0.5)
    lambda_i <- rnorm(N, 0, 0.5)
    f_t <- rnorm(TT, 0, 0.5)

    Y_vec <- numeric(N * TT)
    D_vec <- integer(N * TT)
    id_vec <- integer(N * TT)
    time_vec <- integer(N * TT)

    idx <- 1
    for (i in 1:N) {
        for (t in 1:TT) {
            treated <- (i <= Ntr) && (t > T0)
            D_vec[idx] <- as.integer(treated)
            eps <- rnorm(1, 0, 1)
            Y_vec[idx] <- alpha_i[i] + xi_t[t] +
                lambda_i[i] * f_t[t] +
                tau * D_vec[idx] + eps
            id_vec[idx] <- i
            time_vec[idx] <- t
            idx <- idx + 1
        }
    }

    simdf <- data.frame(
        id = id_vec,
        time = time_vec,
        Y = Y_vec,
        D = D_vec
    )

    out_raw <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "ife",
        r = 1,
        CV = FALSE,
        force = "two-way",
        normalize = FALSE,
        se = FALSE,
        parallel = FALSE
    ))

    out_norm <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "ife",
        r = 1,
        CV = FALSE,
        force = "two-way",
        normalize = TRUE,
        se = FALSE,
        parallel = FALSE
    ))

    ## ATT should be essentially the same
    expect_equal(out_raw$att.avg, out_norm$att.avg, tolerance = 0.05,
        label = "IFE ATT should match with/without normalization")

    ## sigma2.fect should be close
    if (!is.null(out_raw$sigma2.fect) && !is.null(out_norm$sigma2.fect)) {
        rel_diff <- abs(out_raw$sigma2.fect - out_norm$sigma2.fect) /
            max(abs(out_raw$sigma2.fect), 1e-10)
        expect_lt(rel_diff, 0.15,
            label = paste0("IFE sigma2.fect relative difference = ",
                           round(rel_diff, 6)))
    }
})
