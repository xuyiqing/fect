## ============================================================================
## Coverage study for fect()'s built-in ci.method = c("normal", "basic")
## machinery, exercised through fit$est.avg directly (NOT via estimand()).
##
## Same DGPs as run_minimal_coverage.R Scenarios A, B, C; same outer-loop
## parallel design (future workers = 16, K = 200).  For each scenario we fit
## fect() twice --- once with ci.method = "normal" and once with
## ci.method = "basic" --- and read the CI from fit$est.avg.  Coverage =
## fraction of K reps whose CI contains the realized treatment effect (A/B)
## or the population ATT = 3 (C).
##
## Why a separate file: this validates the v2.4.2 ci.method addition on
## fect()'s built-in CI machinery (the path that est.* slots use).  The
## existing run_minimal_coverage.R validates estimand()'s post-hoc CI
## machinery on the same fits.  The two surfaces should agree byte-equally
## on the (normal, basic) overlap; this script is the fect-side proof.
## ============================================================================

suppressPackageStartupMessages({
    library(devtools)
    library(future)
    library(future.apply)
})

.script_path <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0L) {
        return(normalizePath(sub("^--file=", "", file_arg[1])))
    }
    of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(of)) return(normalizePath(of))
    NA_character_
}
.this <- .script_path()
.this_dir <- if (!is.na(.this)) dirname(.this) else getwd()
pkg_root <- normalizePath(file.path(.this_dir, "..", ".."))
if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    stop("Could not locate fect repo root (got '", pkg_root, "').", call. = FALSE)
}
setwd(pkg_root)
cat("Loading fect from:", pkg_root, "\n")
suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))

## DGP helpers (inlined; identical to run_minimal_coverage.R copies so future
## workers don't need to source that file).  Conventions match
## gsynth-note's simulate-xu-{iid,ar1}-rfit2.R: lambda_i ~ U(-sqrt(3), sqrt(3)),
## f_t ~ N(0,1), alpha_i ~ U(-sqrt(3), sqrt(3)), xi_t ~ N(0,1), mu = 5;
## treated cells get ATT_t = t plus N(0, D_sd) heterogeneity.

simulate_factor <- function(seed, N_tr, N_co, TT, T0, r = 2, D_sd = 1, mu = 5,
                            ar_rho = 0.0) {
    set.seed(seed)
    N      <- N_tr + N_co
    T_post <- TT - T0
    ss     <- sqrt(3)

    lambda <- matrix(runif(N * r, -ss, ss), N, r)
    f_mat  <- matrix(rnorm(TT * r), TT, r)
    alpha  <- runif(N, -ss, ss)
    xi     <- rnorm(TT)

    if (ar_rho == 0) {
        e <- matrix(rnorm(TT * N), TT, N)
    } else {
        e <- matrix(NA_real_, TT, N)
        e[1, ] <- rnorm(N)
        shock_sd <- sqrt(1 - ar_rho^2)
        for (t in 2:TT) e[t, ] <- ar_rho * e[t - 1, ] + rnorm(N, 0, shock_sd)
    }

    D <- matrix(0L, TT, N); D[(T0 + 1):TT, 1:N_tr] <- 1L
    eff <- matrix(0, TT, N)
    if (D_sd > 0) {
        eff[(T0 + 1):TT, ] <- matrix(seq_len(T_post), T_post, N) +
                              matrix(rnorm(T_post * N, 0, D_sd), T_post, N)
    } else {
        eff[(T0 + 1):TT, ] <- matrix(seq_len(T_post), T_post, N)
    }

    Y0 <- e + mu + f_mat %*% t(lambda) +
          matrix(alpha, TT, N, byrow = TRUE) +
          matrix(xi,    TT, N, byrow = FALSE)
    Y  <- Y0 + eff * D
    target_att <- mean(eff[(T0 + 1):TT, 1:N_tr])
    panel <- data.frame(
        id = rep(seq_len(N), each = TT),
        time = rep(seq_len(TT), times = N),
        Y = c(Y), D = c(D)
    )
    list(panel = panel, target_att = target_att)
}

simulate_twfe <- function(seed, N_tr, N_co, TT, T0, ar_rho, ATT = 3, mu = 5) {
    set.seed(seed)
    N <- N_tr + N_co
    ss <- sqrt(3)
    alpha <- runif(N, -ss, ss)
    xi    <- rnorm(TT)
    e <- matrix(NA_real_, TT, N)
    e[1, ] <- rnorm(N)
    shock_sd <- sqrt(1 - ar_rho^2)
    for (t in 2:TT) e[t, ] <- ar_rho * e[t - 1, ] + rnorm(N, 0, shock_sd)
    D <- matrix(0L, TT, N); D[(T0 + 1):TT, 1:N_tr] <- 1L
    Y <- e + mu + matrix(alpha, TT, N, byrow = TRUE) +
                  matrix(xi,    TT, N, byrow = FALSE) +
                  ATT * D
    panel <- data.frame(
        id = rep(seq_len(N), each = TT),
        time = rep(seq_len(TT), times = N),
        Y = c(Y), D = c(D)
    )
    list(panel = panel, target_att = ATT)
}

## Per-rep worker: fit fect, extract est.avg row at the requested ci.method
run_one_rep_fect_direct <- function(seed, scenario, ci.method, nboots, pkg_root) {
    ## NOTE: future workers spawn with an installed (CRAN/dev) fect already
    ## resolvable via namespace shipping, so loadedNamespaces() returns
    ## TRUE even on a fresh worker -- but the resolved fect lacks the
    ## v2.4.2 ci.method argument. Force load_all unconditionally on every
    ## call (cheap on warm workers; correct on cold ones).
    if (!isTRUE(getOption(".fect_ci_method_loaded"))) {
        suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
        options(.fect_ci_method_loaded = TRUE)
    }
    t0 <- Sys.time()

    if (scenario == "A") {
        dgp <- simulate_factor(seed, N_tr = 5, N_co = 50, TT = 30, T0 = 20,
                                r = 2, ar_rho = 0)
        fit <- tryCatch(
            suppressMessages(suppressWarnings(
                fect(Y ~ D, data = dgp$panel, index = c("id", "time"),
                     method = "ife", r = 2, force = "two-way", CV = FALSE,
                     se = TRUE, vartype = "parametric", para.error = "auto",
                     ci.method = ci.method,
                     nboots = nboots, parallel = FALSE,
                     time.component.from = "nevertreated")
            )),
            error = function(e) NULL
        )
    } else if (scenario == "B") {
        dgp <- simulate_factor(seed, N_tr = 5, N_co = 50, TT = 30, T0 = 20,
                                r = 2, ar_rho = 0.8)
        fit <- tryCatch(
            suppressMessages(suppressWarnings(
                fect(Y ~ D, data = dgp$panel, index = c("id", "time"),
                     method = "ife", r = 2, force = "two-way", CV = FALSE,
                     se = TRUE, vartype = "parametric", para.error = "auto",
                     ci.method = ci.method,
                     nboots = nboots, parallel = FALSE,
                     time.component.from = "nevertreated")
            )),
            error = function(e) NULL
        )
    } else if (scenario == "C_boot") {
        dgp <- simulate_twfe(seed, N_tr = 20, N_co = 80, TT = 30, T0 = 20,
                              ar_rho = 0.5)
        fit <- tryCatch(
            suppressMessages(suppressWarnings(
                fect(Y ~ D, data = dgp$panel, index = c("id", "time"),
                     method = "fe", force = "two-way", CV = FALSE,
                     se = TRUE, vartype = "bootstrap",
                     ci.method = ci.method,
                     nboots = nboots, parallel = FALSE,
                     time.component.from = "notyettreated")
            )),
            error = function(e) NULL
        )
    } else {
        stop("Unknown scenario: ", scenario)
    }

    fit_wall <- as.numeric(Sys.time() - t0, units = "secs")

    if (is.null(fit) || is.null(fit$est.avg)) {
        return(data.frame(seed = seed, scenario = scenario,
                          ci.method = ci.method,
                          att_hat = NA_real_, se = NA_real_,
                          ci_lo = NA_real_, ci_hi = NA_real_,
                          target_att = dgp$target_att,
                          cover = NA_integer_, width = NA_real_,
                          wall_sec = fit_wall, stringsAsFactors = FALSE))
    }
    row <- fit$est.avg[1, , drop = TRUE]
    att_hat <- unname(row["ATT.avg"])
    ci_lo   <- unname(row["CI.lower"])
    ci_hi   <- unname(row["CI.upper"])
    se      <- unname(row["S.E."])
    cover   <- as.integer(ci_lo <= dgp$target_att && dgp$target_att <= ci_hi)
    data.frame(
        seed       = seed,
        scenario   = scenario,
        ci.method  = ci.method,
        att_hat    = att_hat,
        se         = se,
        ci_lo      = ci_lo,
        ci_hi      = ci_hi,
        target_att = dgp$target_att,
        cover      = cover,
        width      = ci_hi - ci_lo,
        wall_sec   = fit_wall,
        stringsAsFactors = FALSE
    )
}

run_scenario_fect <- function(scenario, label, K, nboots, workers, base_seed,
                               ci_methods = c("normal", "basic")) {
    rows <- list()
    cat(sprintf("\n=== %s ===\n", label))
    cat(sprintf("K=%d reps, nboots=%d, workers=%d (outer parallel), ci.methods=%s\n",
                K, nboots, workers, paste(ci_methods, collapse = ",")))
    seeds <- base_seed + seq_len(K)
    t0 <- Sys.time()
    for (m in ci_methods) {
        results <- future.apply::future_lapply(
            seeds,
            run_one_rep_fect_direct,
            scenario = scenario, ci.method = m,
            nboots = nboots, pkg_root = pkg_root,
            future.seed = TRUE
        )
        rows[[m]] <- do.call(rbind, results)
    }
    df <- do.call(rbind, rows)
    elapsed <- as.numeric(Sys.time() - t0, units = "mins")

    summary_df <- do.call(rbind, lapply(ci_methods, function(m) {
        d <- df[df$ci.method == m, ]
        ok <- !is.na(d$cover)
        n <- sum(ok)
        if (n == 0L) return(NULL)
        cov <- mean(d$cover[ok])
        data.frame(
            scenario     = scenario,
            ci.method    = m,
            n_reps       = n,
            coverage     = cov,
            mc_se        = sqrt(cov * (1 - cov) / n),
            mean_width   = mean(d$width[ok]),
            mean_se      = mean(d$se[ok]),
            empirical_sd = sd(d$att_hat[ok] - d$target_att[ok]),
            mean_bias    = mean(d$att_hat[ok] - d$target_att[ok]),
            wall_min     = round(elapsed, 3),
            stringsAsFactors = FALSE
        )
    }))
    cat(sprintf("  wall: %.2f min\n", elapsed))
    print(summary_df[ , c("ci.method", "n_reps", "coverage", "mc_se",
                          "mean_width", "mean_se", "empirical_sd",
                          "mean_bias")],
          row.names = FALSE, digits = 4)

    list(scenario = scenario, label = label, df = df,
         summary = summary_df, wall_min = elapsed)
}

run_all_fect <- function(K = 200, nboots = 200, workers = 16,
                          out_dir = "/tmp/fect-coverage-study") {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
    future::plan(future::multisession, workers = workers)
    options(future.globals.maxSize = 2e9)
    on.exit(future::plan(future::sequential), add = TRUE)

    ## ci.method in {"normal", "basic"} for all three scenarios.  The
    ## v2.4.2 location-shift fix in R/boot.R (around line 3590) makes
    ## ci.method = "basic" work on parametric fits at the avg-level + per-
    ## event-time CI sites, byte-equivalent to estimand(fit, "att", "basic").
    sA <- run_scenario_fect("A",      "A: factor IID parametric (fect direct)",
                             K, nboots, workers, base_seed = 1000L,
                             ci_methods = c("normal", "basic"))
    sB <- run_scenario_fect("B",      "B: factor AR(1) rho=0.8 parametric (fect direct)",
                             K, nboots, workers, base_seed = 2000L,
                             ci_methods = c("normal", "basic"))
    sC1 <- run_scenario_fect("C_boot", "C1: TWFE AR(1) rho=0.5 bootstrap (fect direct)",
                              K, nboots, workers, base_seed = 3000L,
                              ci_methods = c("normal", "basic"))

    all <- list(A = sA, B = sB, C_boot = sC1)
    for (nm in names(all)) {
        csv <- file.path(out_dir,
                          sprintf("fect_ci_method_%s_K%d_nb%d_%s.csv",
                                  nm, K, nboots, ts))
        write.csv(all[[nm]]$df, csv, row.names = FALSE)
        cat(sprintf("  wrote %s\n", csv))
    }
    summary_df <- do.call(rbind, lapply(all, `[[`, "summary"))
    sum_csv <- file.path(out_dir,
                         sprintf("fect_ci_method_summary_K%d_nb%d_%s.csv",
                                 K, nboots, ts))
    write.csv(summary_df, sum_csv, row.names = FALSE)

    cat("\n", strrep("=", 78), "\n", sep = "")
    cat("FECT CI.METHOD COVERAGE SUMMARY (K=", K, ", nboots=", nboots, ")\n",
        sep = "")
    cat(strrep("=", 78), "\n", sep = "")
    print(summary_df[, c("scenario", "ci.method", "coverage", "mc_se",
                         "mean_width", "mean_se", "empirical_sd")],
          row.names = FALSE, digits = 4)
    cat(sprintf("\nSummary CSV: %s\n", sum_csv))
    cat(strrep("=", 78), "\n", sep = "")

    invisible(list(scenarios = all, summary = summary_df,
                   summary_csv = sum_csv))
}

if (sys.nframe() == 0L || identical(environmentName(parent.frame()), "R_GlobalEnv")) {
    run_all_fect(K = 200, nboots = 200, workers = 16)
}
