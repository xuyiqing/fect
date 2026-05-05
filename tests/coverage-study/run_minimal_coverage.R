## ============================================================================
## Minimal coverage validation for the v2.4.2 inference paths.
##
## Three scenarios, no covariates throughout:
##
##   A. Factor model (r=2), IID errors, parametric inference
##   B. Factor model (r=2), AR(1) rho=0.8 errors, parametric inference
##   C. Large-N additive TWFE (r=0), AR(1) rho=0.5 errors,
##      bootstrap (cluster) AND jackknife inference (two cells, same DGP)
##
## DGPs A and B replicate gsynth-note's canonical Xu-2017 setup
## (code/sims/coverage/simulate-xu-{iid,ar1}-rfit2.R). DGP C scales N
## up to test the cluster-bootstrap and jackknife paths under serial
## correlation that fect's empirical-residual draw can capture.
##
## Parallelization: outer-loop (across reps) via future + future.apply.
## Each fect() call runs sequentially (parallel = FALSE) inside a worker.
## Worker pool = 16; reps are embarrassingly parallel so this beats
## inner bootstrap parallelism (no per-rep cluster fork/join overhead).
##
## Wall time: ~5-8 min on 16 cores.  Output:
##   /tmp/fect-coverage-study/minimal_<scenario>_K<K>_nb<nb>_<ts>.csv
## ============================================================================

suppressPackageStartupMessages({
    library(devtools)
    library(future)
    library(future.apply)
})

## resolve repo root robustly (Rscript: --file=<path>; source(): ofile)
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
.this_file <- .script_path()
script_dir <- if (!is.na(.this_file)) dirname(.this_file) else getwd()
pkg_root   <- normalizePath(file.path(script_dir, "..", ".."))
if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    stop("Could not locate fect repo root (got '", pkg_root,
         "'); expected DESCRIPTION there. Run from repo root or pass ",
         "--file=<absolute path>.", call. = FALSE)
}
setwd(pkg_root)
cat("Loading fect from:", pkg_root, "\n")
suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))

## ---- DGP helpers (identical structure across scenarios) -------------------
##
## Conventions match gsynth-note simulate-xu-*.R:
##   loadings lambda_i ~ U(-sqrt(3), sqrt(3))   variance 1
##   factors  f_t      ~ N(0, 1)
##   unit FE  alpha_i  ~ U(-sqrt(3), sqrt(3))   variance 1
##   time FE  xi_t     ~ N(0, 1)
##   grand mean mu = 5
##   Treatment block fixed at units 1:N_tr, periods (T0+1):TT
##   ATT_t  = t for t in 1:T_post, plus N(0, D_sd) per (unit, time) in treated post
##   Coverage target = realized average treated-post effect (within rep).

simulate_factor <- function(seed, N_tr, N_co, TT, T0, r = 2, D_sd = 1, mu = 5,
                            ar_rho = 0.0) {
    set.seed(seed)
    N      <- N_tr + N_co
    T_post <- TT - T0
    ss     <- sqrt(3)

    lambda <- matrix(runif(N * r, -ss, ss), N, r)
    f_mat  <- matrix(rnorm(TT * r),         TT, r)
    alpha  <- runif(N,  -ss, ss)
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
        id   = rep(seq_len(N), each  = TT),
        time = rep(seq_len(TT), times = N),
        Y    = c(Y),
        D    = c(D)
    )
    list(panel = panel, target_att = target_att)
}

## Additive TWFE (r=0), constant ATT.  Used by Scenario C.
simulate_twfe <- function(seed, N_tr, N_co, TT, T0, ar_rho, ATT = 3, mu = 5) {
    set.seed(seed)
    N <- N_tr + N_co
    ss <- sqrt(3)

    alpha <- runif(N,  -ss, ss)
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
        id   = rep(seq_len(N), each  = TT),
        time = rep(seq_len(TT), times = N),
        Y    = c(Y),
        D    = c(D)
    )
    list(panel = panel, target_att = ATT)
}

## ---- Per-scenario fit + extract ATT/CI ------------------------------------
## Each fit_* function returns the fect fit object (or NULL on error).
## All fect() calls run with parallel = FALSE (outer parallelism handles K).
## keep.sims = TRUE on parametric / bootstrap so estimand() can dispatch
## across all 5 ci.methods on the same fit (one fit per rep, not five).

fit_parametric_ife <- function(panel, nboots) {
    tryCatch(
        suppressMessages(suppressWarnings(
            fect(Y ~ D, data = panel, index = c("id", "time"),
                 method = "ife", r = 2, force = "two-way",
                 CV = FALSE,
                 se = TRUE, vartype = "parametric", para.error = "auto",
                 nboots = nboots, parallel = FALSE,
                 keep.sims = TRUE,
                 time.component.from = "nevertreated")
        )),
        error = function(e) NULL
    )
}

fit_bootstrap_fe <- function(panel, nboots) {
    tryCatch(
        suppressMessages(suppressWarnings(
            fect(Y ~ D, data = panel, index = c("id", "time"),
                 method = "fe", force = "two-way",
                 CV = FALSE,
                 se = TRUE, vartype = "bootstrap",
                 nboots = nboots, parallel = FALSE,
                 keep.sims = TRUE,
                 time.component.from = "notyettreated")
        )),
        error = function(e) NULL
    )
}

fit_jackknife_fe <- function(panel) {
    tryCatch(
        suppressMessages(suppressWarnings(
            fect(Y ~ D, data = panel, index = c("id", "time"),
                 method = "fe", force = "two-way",
                 CV = FALSE,
                 se = TRUE, vartype = "jackknife",
                 parallel = FALSE,
                 time.component.from = "notyettreated")
        )),
        error = function(e) NULL
    )
}

## Pull (estimate, ci.lo, ci.hi, se) from a fit at a specific ci.method.
## Uses estimand(fit, "att", "overall") to dispatch by ci.method on
## parametric and bootstrap fits; reads fit$est.avg directly for jackknife
## (estimand path requires keep.sims = TRUE which jackknife does not yield).
extract_at <- function(fit, ci.method) {
    null_row <- list(att_hat = NA_real_, se = NA_real_,
                     ci_lo = NA_real_, ci_hi = NA_real_, ok = FALSE)
    if (is.null(fit)) return(null_row)
    if (isTRUE(fit$vartype == "jackknife")) {
        if (is.null(fit$est.avg)) return(null_row)
        row <- fit$est.avg[1, , drop = TRUE]
        return(list(att_hat = unname(row["ATT.avg"]),
                    se      = unname(row["S.E."]),
                    ci_lo   = unname(row["CI.lower"]),
                    ci_hi   = unname(row["CI.upper"]),
                    ok      = TRUE))
    }
    est <- tryCatch(
        suppressMessages(suppressWarnings(
            estimand(fit, type = "att", by = "overall",
                     ci.method = ci.method)
        )),
        error = function(e) NULL
    )
    if (is.null(est) || nrow(est) == 0L) return(null_row)
    row <- est[1, , drop = TRUE]
    list(att_hat = unname(row[["estimate"]]),
         se      = unname(row[["se"]]),
         ci_lo   = unname(row[["ci.lo"]]),
         ci_hi   = unname(row[["ci.hi"]]),
         ok      = TRUE)
}

## ---- Worker function ------------------------------------------------------
## Runs once per (seed, scenario).  First call in each worker triggers
## load_all(); subsequent calls reuse the loaded namespace.
## Returns one row per (rep, ci.method).  Jackknife scenario forces
## ci.methods = "normal" (only valid CI for jackknife per E&T 1993 Ch 11).

CI_METHODS_FULL <- c("basic", "percentile", "bc", "bca", "normal")

run_one_rep <- function(seed, scenario, nboots, pkg_root,
                         ci_methods = CI_METHODS_FULL) {
    if (!"fect" %in% loadedNamespaces()) {
        suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
    }
    t0 <- Sys.time()

    if (scenario == "A") {
        dgp <- simulate_factor(seed, N_tr = 5, N_co = 50, TT = 30, T0 = 20,
                                r = 2, ar_rho = 0)
        fit <- fit_parametric_ife(dgp$panel, nboots)
        cms <- ci_methods
    } else if (scenario == "B") {
        dgp <- simulate_factor(seed, N_tr = 5, N_co = 50, TT = 30, T0 = 20,
                                r = 2, ar_rho = 0.8)
        fit <- fit_parametric_ife(dgp$panel, nboots)
        cms <- ci_methods
    } else if (scenario == "C_boot") {
        dgp <- simulate_twfe(seed, N_tr = 20, N_co = 80, TT = 30, T0 = 20,
                              ar_rho = 0.5)
        fit <- fit_bootstrap_fe(dgp$panel, nboots)
        cms <- ci_methods
    } else if (scenario == "C_jack") {
        dgp <- simulate_twfe(seed, N_tr = 20, N_co = 80, TT = 30, T0 = 20,
                              ar_rho = 0.5)
        fit <- fit_jackknife_fe(dgp$panel)
        cms <- "normal"          # only valid CI for jackknife (E&T 1993 §11)
    } else {
        stop("Unknown scenario: ", scenario)
    }

    fit_wall <- as.numeric(Sys.time() - t0, units = "secs")

    rows <- lapply(cms, function(m) {
        out <- extract_at(fit, m)
        cover <- if (out$ok) {
            as.integer(out$ci_lo <= dgp$target_att &&
                       dgp$target_att <= out$ci_hi)
        } else NA_integer_
        data.frame(
            seed       = seed,
            scenario   = scenario,
            ci.method  = m,
            att_hat    = out$att_hat,
            se         = out$se,
            ci_lo      = out$ci_lo,
            ci_hi      = out$ci_hi,
            target_att = dgp$target_att,
            cover      = cover,
            width      = out$ci_hi - out$ci_lo,
            wall_sec   = fit_wall,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

## ---- Driver: parallelize outer loop across reps ---------------------------

summarise_per_method <- function(df) {
    parts <- split(df, df$ci.method)
    out <- lapply(names(parts), function(m) {
        d  <- parts[[m]]
        ok <- !is.na(d$cover)
        n  <- sum(ok)
        if (n == 0L) return(NULL)
        cov     <- mean(d$cover[ok])
        mc_se   <- sqrt(cov * (1 - cov) / n)
        data.frame(
            ci.method    = m,
            n_reps       = n,
            coverage     = cov,
            mc_se        = mc_se,
            mean_width   = mean(d$width[ok]),
            mean_se      = mean(d$se[ok], na.rm = TRUE),
            empirical_sd = sd(d$att_hat[ok] - d$target_att[ok]),
            mean_bias    = mean(d$att_hat[ok] - d$target_att[ok]),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out[!sapply(out, is.null)])
}

run_scenario <- function(scenario, label, K, nboots, workers, base_seed,
                          ci_methods = CI_METHODS_FULL) {
    cat(sprintf("\n=== %s ===\n", label))
    cat(sprintf("K=%d reps, nboots=%d, workers=%d, ci.methods=%s\n",
                K, nboots, workers, paste(ci_methods, collapse = ",")))

    seeds <- base_seed + seq_len(K)
    t0 <- Sys.time()

    results <- future.apply::future_lapply(
        seeds,
        run_one_rep,
        scenario   = scenario,
        nboots     = nboots,
        pkg_root   = pkg_root,
        ci_methods = ci_methods,
        future.seed = TRUE
    )
    df <- do.call(rbind, results)
    elapsed <- as.numeric(Sys.time() - t0, units = "mins")

    summary_df <- summarise_per_method(df)
    summary_df <- cbind(scenario = scenario, summary_df,
                        wall_min = round(elapsed, 3))

    cat(sprintf("  wall: %.2f min\n", elapsed))
    print(summary_df[ , c("ci.method", "n_reps", "coverage", "mc_se",
                          "mean_width", "mean_se", "empirical_sd",
                          "mean_bias")],
          row.names = FALSE, digits = 4)

    list(scenario = scenario, label = label, df = df,
         summary = summary_df, wall_min = elapsed)
}

run_all <- function(K = 200, nboots = 200, workers = 16,
                    ci_methods = CI_METHODS_FULL,
                    out_dir = "/tmp/fect-coverage-study") {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ts <- format(Sys.time(), "%Y%m%d-%H%M%S")

    future::plan(future::multisession, workers = workers)
    options(future.globals.maxSize = 2e9)
    on.exit(future::plan(future::sequential), add = TRUE)

    sA <- run_scenario("A",      "A: factor IID parametric",
                       K, nboots, workers, base_seed = 1000L,
                       ci_methods = ci_methods)
    sB <- run_scenario("B",      "B: factor AR(1) rho=0.8 parametric",
                       K, nboots, workers, base_seed = 2000L,
                       ci_methods = ci_methods)
    sC1 <- run_scenario("C_boot", "C1: TWFE AR(1) rho=0.5 bootstrap",
                        K, nboots, workers, base_seed = 3000L,
                        ci_methods = ci_methods)
    sC2 <- run_scenario("C_jack", "C2: TWFE AR(1) rho=0.5 jackknife",
                        K, nboots, workers, base_seed = 4000L,
                        ci_methods = "normal")

    all <- list(A = sA, B = sB, C_boot = sC1, C_jack = sC2)
    for (nm in names(all)) {
        csv <- file.path(out_dir,
                          sprintf("minimal_%s_K%d_nb%d_%s.csv",
                                  nm, K, nboots, ts))
        write.csv(all[[nm]]$df, csv, row.names = FALSE)
        cat(sprintf("  wrote %s\n", csv))
    }

    summary_df <- do.call(rbind, lapply(all, `[[`, "summary"))
    sum_csv <- file.path(out_dir,
                         sprintf("minimal_summary_K%d_nb%d_%s.csv",
                                 K, nboots, ts))
    write.csv(summary_df, sum_csv, row.names = FALSE)

    cat("\n", strrep("=", 78), "\n", sep = "")
    cat("MINIMAL COVERAGE SUMMARY  (K=", K, ", nboots=", nboots, ")\n",
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

## Run when invoked as a top-level Rscript (not when sourced for definitions
## by run_minimal_coverage_tail_rerun.R; that script sets the option below).
if (!isTRUE(getOption(".fect_minimal_no_autorun")) &&
    (sys.nframe() == 0L ||
     identical(environmentName(parent.frame()), "R_GlobalEnv"))) {
    run_all(K = 200, nboots = 200, workers = 16)
}
