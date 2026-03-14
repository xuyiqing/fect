## ---------------------------------------------------------
## Simulation test: Staggered adoption design
## DGP: Three treatment cohorts with tau=2.0
## ---------------------------------------------------------

test_that("FE estimator recovers ATT under staggered adoption", {
    skip_on_cran()

    set.seed(5001)

    N <- 30
    TT <- 20
    tau <- 2.0
    nsim <- 50

    ## Three cohorts: treated at t=8, t=12, t=16
    ## 10 units each: first 10 = cohort 1, next 10 = cohort 2, last 10 = control
    cohort_times <- c(rep(8, 10), rep(12, 10), rep(Inf, 10))

    estimates <- numeric(nsim)

    for (sim in 1:nsim) {
        alpha_i <- rnorm(N, 0, 1)
        xi_t <- rnorm(TT, 0, 0.5)

        Y_vec <- numeric(N * TT)
        D_vec <- integer(N * TT)
        id_vec <- integer(N * TT)
        time_vec <- integer(N * TT)

        idx <- 1
        for (i in 1:N) {
            for (t in 1:TT) {
                treated <- (t >= cohort_times[i])
                D_vec[idx] <- as.integer(treated)
                eps <- rnorm(1, 0, 1)
                Y_vec[idx] <- alpha_i[i] + xi_t[t] +
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
        label = paste0("Staggered FE bias = ", round(bias, 4),
                       ", mean_est = ", round(mean_est, 4),
                       ", tau = ", tau))

    sd_est <- sd(estimates, na.rm = TRUE)
    expect_lt(sd_est, 1.5,
        label = paste0("Staggered FE SD = ", round(sd_est, 4)))
})

test_that("IFE handles staggered adoption with factor structure", {
    skip_on_cran()

    set.seed(5002)

    N <- 30
    TT <- 20
    tau <- 2.0

    ## Three cohorts: treated at t=8, t=12, t=16
    cohort_times <- c(rep(8, 10), rep(12, 10), rep(Inf, 10))

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
            treated <- (t >= cohort_times[i])
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

    ## Check that ATT is in a reasonable range
    expect_lt(abs(out$att.avg - tau), 1.5,
        label = paste0("Staggered IFE: att.avg = ",
                       round(out$att.avg, 4), ", tau = ", tau))

    ## Pre-treatment effects should be small
    pre_idx <- which(out$time <= 0)
    if (length(pre_idx) > 0) {
        pre_effects <- out$att[pre_idx]
        pre_effects <- pre_effects[!is.na(pre_effects)]
        if (length(pre_effects) > 0) {
            mean_pre <- mean(abs(pre_effects))
            expect_lt(mean_pre, 2.0,
                label = paste0("Staggered IFE |pre-treatment effect| = ",
                               round(mean_pre, 4)))
        }
    }
})
