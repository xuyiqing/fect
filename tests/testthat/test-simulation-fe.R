## ---------------------------------------------------------
## Simulation test: Fixed Effects (FE) estimator
## DGP: Two-way FE model with known treatment effect tau=3.0
## ---------------------------------------------------------

test_that("FE estimator recovers ATT under simple two-way FE DGP", {

  skip_on_cran()

    set.seed(2024)

    N <- 30       # units
    TT <- 15      # time periods
    T0 <- 10      # treatment starts at period 11 for treated units
    Ntr <- 10     # number of treated units
    tau <- 3.0    # true treatment effect
    nsim <- 50    # number of Monte Carlo replications

    estimates <- numeric(nsim)

    for (sim in 1:nsim) {
        ## Unit and time fixed effects
        alpha_i <- rnorm(N, 0, 1)
        xi_t <- rnorm(TT, 0, 0.5)

        ## Generate panel data
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

        out <- suppressWarnings(fect::fect(
            Y ~ D,
            data = simdf,
            index = c("id", "time"),
            method = "fe",
            force = "two-way",
            se = FALSE,
            parallel = FALSE
        ))

        estimates[sim] <- out$att.avg
    }

    ## Check bias: mean estimate should be close to true tau
    mean_est <- mean(estimates, na.rm = TRUE)
    bias <- abs(mean_est - tau)
    expect_lt(bias, 0.5,
        label = paste0("FE bias = ", round(bias, 4),
                       ", mean_est = ", round(mean_est, 4),
                       ", tau = ", tau))

    ## Check dispersion: SD should be reasonable
    sd_est <- sd(estimates, na.rm = TRUE)
    expect_lt(sd_est, 1.5,
        label = paste0("FE SD = ", round(sd_est, 4)))

    ## Check that estimates are centered around tau
    ## (median should also be near tau)
    median_est <- median(estimates, na.rm = TRUE)
    expect_lt(abs(median_est - tau), 0.6,
        label = paste0("FE median bias = ", round(abs(median_est - tau), 4)))
})

test_that("FE estimator produces near-zero pre-treatment effects", {

  skip_on_cran()

    set.seed(2025)

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

    out <- suppressWarnings(fect::fect(
        Y ~ D,
        data = simdf,
        index = c("id", "time"),
        method = "fe",
        force = "two-way",
        se = FALSE,
        parallel = FALSE
    ))

    ## Pre-treatment dynamic effects should be near zero
    pre_idx <- which(out$time <= 0)
    if (length(pre_idx) > 0) {
        pre_effects <- out$att[pre_idx]
        pre_effects <- pre_effects[!is.na(pre_effects)]
        if (length(pre_effects) > 0) {
            mean_pre <- mean(abs(pre_effects))
            expect_lt(mean_pre, 1.5,
                label = paste0("Mean |pre-treatment effect| = ",
                               round(mean_pre, 4)))
        }
    }
})
