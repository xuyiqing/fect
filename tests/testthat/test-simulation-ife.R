## ---------------------------------------------------------
## Simulation test: Interactive Fixed Effects (IFE) estimator
## DGP: Factor model with one factor and tau=2.0
## ---------------------------------------------------------

test_that("IFE estimator recovers ATT under factor model DGP", {

  skip_on_cran()

    set.seed(3001)

    N <- 30
    TT <- 15
    T0 <- 10
    Ntr <- 10
    tau <- 2.0
    nsim <- 50

    estimates <- numeric(nsim)

    for (sim in 1:nsim) {
        ## Unit and time fixed effects
        alpha_i <- rnorm(N, 0, 1)
        xi_t <- rnorm(TT, 0, 0.5)

        ## Single factor structure
        lambda_i <- rnorm(N, 0, 1)
        f_t <- rnorm(TT, 0, 1)

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

        out <- suppressWarnings(fect::fect(
            Y ~ D,
            data = simdf,
            index = c("id", "time"),
            method = "ife",
            r = 1,
            CV = FALSE,
            force = "two-way",
            se = FALSE,
            parallel = FALSE
        ))

        estimates[sim] <- out$att.avg
    }

    mean_est <- mean(estimates, na.rm = TRUE)
    bias <- abs(mean_est - tau)
    expect_lt(bias, 0.8,
        label = paste0("IFE bias = ", round(bias, 4),
                       ", mean_est = ", round(mean_est, 4),
                       ", tau = ", tau))

    sd_est <- sd(estimates, na.rm = TRUE)
    expect_lt(sd_est, 2.0,
        label = paste0("IFE SD = ", round(sd_est, 4)))
})

test_that("IFE with r=0 reduces to FE and produces consistent estimates", {

  skip_on_cran()

    set.seed(3002)

    N <- 30
    TT <- 15
    T0 <- 10
    Ntr <- 10
    tau <- 2.5

    ## Pure two-way FE DGP (no factors)
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
            eps <- rnorm(1, 0, 0.5)
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

    out_fe <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "fe",
        force = "two-way",
        se = FALSE,
        parallel = FALSE
    ))

    out_ife <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "ife",
        r = 0,
        CV = FALSE,
        force = "two-way",
        se = FALSE,
        parallel = FALSE
    ))

    ## With r=0, IFE should produce the same result as FE
    expect_equal(out_fe$att.avg, out_ife$att.avg, tolerance = 0.01,
        label = "IFE(r=0) should match FE")
})
