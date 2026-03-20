## ---------------------------------------------------------
## Simulation test: Complex Fixed Effects (CFE) estimator
## DGP: Two-way FE with a time-varying covariate and tau=2.5
## ---------------------------------------------------------

test_that("CFE estimator recovers ATT with time-varying covariate", {

  skip_on_cran()

    set.seed(4001)

    N <- 30
    TT <- 15
    T0 <- 10
    Ntr <- 10
    tau <- 2.5
    beta_x <- 1.5  # true coefficient on X
    nsim <- 50

    estimates <- numeric(nsim)

    for (sim in 1:nsim) {
        alpha_i <- rnorm(N, 0, 1)
        xi_t <- rnorm(TT, 0, 0.5)

        Y_vec <- numeric(N * TT)
        D_vec <- integer(N * TT)
        X1_vec <- numeric(N * TT)
        id_vec <- integer(N * TT)
        time_vec <- integer(N * TT)

        idx <- 1
        for (i in 1:N) {
            for (t in 1:TT) {
                treated <- (i <= Ntr) && (t > T0)
                D_vec[idx] <- as.integer(treated)
                X1_vec[idx] <- rnorm(1, 0, 1)
                eps <- rnorm(1, 0, 1)
                Y_vec[idx] <- alpha_i[i] + xi_t[t] +
                    beta_x * X1_vec[idx] +
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
            D = D_vec,
            X1 = X1_vec
        )

        out <- suppressWarnings(fect::fect(
            Y ~ D + X1,
            data = simdf,
            index = c("id", "time"),
            method = "fe",
            force = "two-way",
            se = FALSE,
            parallel = FALSE
        ))

        estimates[sim] <- out$att.avg
    }

    mean_est <- mean(estimates, na.rm = TRUE)
    bias <- abs(mean_est - tau)
    expect_lt(bias, 0.5,
        label = paste0("CFE bias = ", round(bias, 4),
                       ", mean_est = ", round(mean_est, 4),
                       ", tau = ", tau))

    sd_est <- sd(estimates, na.rm = TRUE)
    expect_lt(sd_est, 1.5,
        label = paste0("CFE SD = ", round(sd_est, 4)))
})

test_that("CFE estimator handles correlated covariate correctly", {

  skip_on_cran()

    set.seed(4002)

    N <- 40
    TT <- 15
    T0 <- 10
    Ntr <- 15
    tau <- 2.5
    beta_x <- 2.0

    alpha_i <- rnorm(N, 0, 1)
    xi_t <- rnorm(TT, 0, 0.5)

    Y_vec <- numeric(N * TT)
    D_vec <- integer(N * TT)
    X1_vec <- numeric(N * TT)
    id_vec <- integer(N * TT)
    time_vec <- integer(N * TT)

    idx <- 1
    for (i in 1:N) {
        for (t in 1:TT) {
            treated <- (i <= Ntr) && (t > T0)
            D_vec[idx] <- as.integer(treated)
            ## Covariate correlated with unit FE
            X1_vec[idx] <- 0.5 * alpha_i[i] + rnorm(1, 0, 1)
            eps <- rnorm(1, 0, 1)
            Y_vec[idx] <- alpha_i[i] + xi_t[t] +
                beta_x * X1_vec[idx] +
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
        D = D_vec,
        X1 = X1_vec
    )

    out <- suppressWarnings(fect::fect(
        Y ~ D + X1,
        data = simdf,
        index = c("id", "time"),
        method = "fe",
        force = "two-way",
        se = FALSE,
        parallel = FALSE
    ))

    ## ATT should still be close to true tau
    expect_lt(abs(out$att.avg - tau), 1.0,
        label = paste0("CFE with correlated X: att.avg = ",
                       round(out$att.avg, 4), ", tau = ", tau))

    ## Output structure checks
    expect_true(is.numeric(out$att.avg))
    expect_true(is.matrix(out$eff))
    expect_true(!is.null(out$beta))
})
