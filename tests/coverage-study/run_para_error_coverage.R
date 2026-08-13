## ============================================================================
## Coverage validation suite for the parametric / para.error inference path.
##
## RUN WHEN: any code change touches the bootstrap distribution.  Specifically:
##   - R/boot.R parametric branch (line ~774, the method %in% {gsynth,ife,cfe}
##     parametric path) or its para.error sub-dispatch
##   - R/po-estimands.R location-shift code (search "is_parametric")
##   - the vartype, ci.method, or para.error enum
##   - jackknife dispatch / slot contract
##   - eff.boot or att.avg.boot population
##
## DO NOT RUN ON: doc edits, vignette changes, R/plot.R, refactors that don't
## touch the bootstrap distribution.  The fast `devtools::test()` suite is the
## right gate for those.
##
## Wall time: ~30 min in parallel (T19 + T20 each ~30 min, T21 ~10 min;
## scripts run sequentially below by default for reproducibility).
##
## Output: /tmp/fect-coverage-study/coverage_para_error_<timestamp>.csv
##         + plain-text PASS/FAIL summary printed to stdout.
##
## Acceptance:
##   T19 (DGP-A, IID Gaussian, n=100, B=1000)
##     coverage in [0.90, 0.99] for all (para.error mode x ci.method) cells
##     (0.90 not 0.95: small-N IID parametric bootstrap targets V_t alone
##     and misses Var_{Lambda,F}[b_t]; SE / empirical-SD ~ 0.91 -> ~0.91
##     coverage analytically. See README and 2026-05-03 run log.)
##   T20 (DGP-A8, AR(1) rho=0.8, n=100, B=1000)
##     coverage >= 0.91 for all cells (AR(1) inflates variance)
##   T21 (DGP-A, n=50, B=500)
##     wild/empirical CI width ratio in [0.70, 1.30] across all 5 ci.methods
##
## Background (gsynth-note section 2): the parametric pseudo-treated bootstrap
## targets the conditional variance V_t = Var(ATT_hat - ATT | Lambda, F, X, D).
## Treatment block D is held fixed across replications (D[(T0+1):TT, 1:Ntr] <- 1)
## so the empirical variance across reps tracks E_{(Lambda,F)}[V_t], not the
## marginal variance Var(ATT_hat - ATT) which by the law of total variance equals
## E[V_t] + Var_{(Lambda,F)}[b_t] and includes a finite-sample-bias contribution
## the bootstrap is silent about.  Re-randomizing D would inject Var_D[b_t] and
## bias measured coverage downward.
## ============================================================================

suppressPackageStartupMessages({ library(devtools) })

## resolve repo root from script location (works when called from anywhere)
script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
pkg_root <- normalizePath(file.path(script_dir, "..", ".."))
setwd(pkg_root)
cat("Loading fect from:", pkg_root, "\n")
devtools::load_all(pkg_root, quiet = TRUE)

## ------------------------------------------------------------ DGP helpers --

## DGP-A: additive TWFE, IID Gaussian, true ATT = 3.0, fully observed
## Treatment block fixed at units 1:Ntr (see gsynth-note section 2).
dgp_a <- function(seed) {
    set.seed(seed)
    N <- 40; TT <- 20; T0 <- 12; Ntr <- 12
    alpha_i <- rnorm(N, 0, 2)
    xi_t    <- rnorm(TT, 0, 1)
    eps     <- matrix(rnorm(N * TT, 0, 1), TT, N)
    D <- matrix(0L, TT, N); D[(T0 + 1):TT, 1:Ntr] <- 1L
    Y <- outer(xi_t, rep(1, N)) + outer(rep(1, TT), alpha_i) + 3 * D + eps
    data.frame(id = rep(1:N, each = TT), time = rep(1:TT, N),
               Y  = c(Y), D = c(D))
}

## DGP-A8: same as DGP-A but errors are AR(1) with rho = 0.8
dgp_a8 <- function(seed) {
    set.seed(seed)
    N <- 40; TT <- 20; T0 <- 12; Ntr <- 12
    alpha_i <- rnorm(N, 0, 2)
    xi_t    <- rnorm(TT, 0, 1)
    eps <- matrix(NA, TT, N)
    for (i in 1:N) {
        e <- rnorm(TT, 0, 1)
        eps[1, i] <- e[1] / sqrt(1 - 0.64)   # stationary init
        for (t in 2:TT) eps[t, i] <- 0.8 * eps[t - 1, i] + e[t]
    }
    D <- matrix(0L, TT, N); D[(T0 + 1):TT, 1:Ntr] <- 1L
    Y <- outer(xi_t, rep(1, N)) + outer(rep(1, TT), alpha_i) + 3 * D + eps
    data.frame(id = rep(1:N, each = TT), time = rep(1:TT, N),
               Y  = c(Y), D = c(D))
}

## ------------------------------------------------------------ shared call --

fect_para <- function(df, para.error, nboots = 1000, parallel = TRUE, cores = 10) {
    suppressMessages(suppressWarnings(
        fect(Y ~ D, data = df, index = c("id", "time"),
             method = "fe", force = "two-way",
             se = TRUE, vartype = "parametric", para.error = para.error,
             nboots = nboots, parallel = parallel, cores = cores,
             time.component.from = "nevertreated",
             keep.sims = TRUE, CV = FALSE)
    ))
}

## ------------------------------------------------------------ T19 / T20 ---

run_coverage_dgp <- function(label, dgp_fn, n_reps = 100, nboots = 1000,
                              cores = 10, threshold = 0.91, upper = 0.99,
                              true_att = 3.0) {
    ci_methods <- c("basic", "percentile", "bc", "bca", "normal")
    para_modes <- c("ar", "empirical", "wild")

    raw_covers <- array(NA_real_,
                        dim = c(n_reps, length(para_modes), length(ci_methods)),
                        dimnames = list(NULL, para_modes, ci_methods))
    ci_widths <- raw_covers
    att_hat   <- matrix(NA_real_, n_reps, length(para_modes),
                        dimnames = list(NULL, para_modes))

    cat(sprintf("\n=== %s: %d reps x %d nboots x 3 modes x 5 ci.methods ===\n",
                label, n_reps, nboots))
    t0 <- Sys.time()
    for (r in seq_len(n_reps)) {
        df <- dgp_fn(seed = r * 100)
        for (pm in para_modes) {
            fit <- fect_para(df, para.error = pm, nboots = nboots,
                             parallel = TRUE, cores = cores)
            att_hat[r, pm] <- fit$att.avg
            for (m in ci_methods) {
                est <- estimand(fit, "att", "overall",
                                window = c(1, 8), ci.method = m)
                raw_covers[r, pm, m] <- as.integer(est$ci.lo <= true_att &&
                                                    est$ci.hi >= true_att)
                ci_widths[r, pm, m] <- est$ci.hi - est$ci.lo
            }
        }
        if (r %% 5 == 0) {
            elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
            cat(sprintf("  rep %3d / %d  elapsed %.1f min  eta %.1f min\n",
                        r, n_reps, elapsed, elapsed / r * (n_reps - r)))
        }
    }

    cov_mat   <- apply(raw_covers, c(2, 3), mean, na.rm = TRUE)
    width_mat <- apply(ci_widths,  c(2, 3), mean, na.rm = TRUE)
    bias_vec  <- colMeans(att_hat, na.rm = TRUE) - true_att

    cat(sprintf("\n%s coverage table:\n", label))
    print(round(cov_mat, 3))
    cat("\nMean CI width:\n"); print(round(width_mat, 3))
    cat("\nBias (att.hat - true):\n"); print(round(bias_vec, 3))

    pass <- all(cov_mat >= threshold & cov_mat <= upper)
    cat(sprintf("\n%s PASS at [%.2f, %.2f] threshold: %s\n",
                label, threshold, upper, pass))

    list(label = label, cov_mat = cov_mat, width_mat = width_mat,
         bias = bias_vec, pass = pass, raw_covers = raw_covers,
         ci_widths = ci_widths, att_hat = att_hat,
         n_reps = n_reps, nboots = nboots)
}

## ------------------------------------------------------------ T21 width --

run_width_parity <- function(n_reps = 50, nboots = 500, cores = 10,
                             ratio_lo = 0.70, ratio_hi = 1.30) {
    ci_methods <- c("basic", "percentile", "bc", "bca", "normal")
    modes      <- c("empirical", "wild")
    width <- array(NA_real_, dim = c(n_reps, length(modes), length(ci_methods)),
                   dimnames = list(NULL, modes, ci_methods))

    cat(sprintf("\n=== T21 width parity: %d reps x %d nboots x 2 modes ===\n",
                n_reps, nboots))
    t0 <- Sys.time()
    for (r in seq_len(n_reps)) {
        df <- dgp_a(seed = r * 77)
        for (md in modes) {
            fit <- fect_para(df, para.error = md, nboots = nboots,
                             parallel = TRUE, cores = cores)
            for (m in ci_methods) {
                est <- estimand(fit, "att", "overall",
                                window = c(1, 8), ci.method = m)
                width[r, md, m] <- est$ci.hi - est$ci.lo
            }
        }
        if (r %% 5 == 0) {
            elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
            cat(sprintf("  rep %2d / %d  elapsed %.1f min\n",
                        r, n_reps, elapsed))
        }
    }
    mean_width <- apply(width, c(2, 3), mean, na.rm = TRUE)
    ratios <- mean_width["wild", ] / mean_width["empirical", ]

    cat("\nMean CI width by mode x ci.method:\n"); print(round(mean_width, 3))
    cat("\nratio wild / empirical:\n"); print(round(ratios, 3))
    pass <- all(ratios >= ratio_lo & ratios <= ratio_hi)
    cat(sprintf("\nT21 PASS [%.2f, %.2f]: %s\n", ratio_lo, ratio_hi, pass))
    list(mean_width = mean_width, ratios = ratios, width = width, pass = pass)
}

## ------------------------------------------------------------ entrypoint --

run_all <- function(out_dir = "/tmp/fect-coverage-study") {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ts <- format(Sys.time(), "%Y%m%d-%H%M%S")

    ## T19 threshold = 0.90 (not 0.91): at N=40 IID with no factor structure,
    ## the parametric pseudo-treated bootstrap targets V_t alone and misses
    ## Var_{Lambda,F}[b_t]; bootstrap SE / empirical SD ~ 0.91, yielding ~0.91
    ## coverage analytically (see README "Why fixed treated block" + 2026-05-03
    ## coverage-completion run log). T20 keeps 0.91 (AR(1) inflates variance,
    ## empirically lands at 0.96+).
    t19 <- run_coverage_dgp("T19 (DGP-A, IID)",     dgp_a,  n_reps = 100, nboots = 1000,
                             threshold = 0.90)
    t20 <- run_coverage_dgp("T20 (DGP-A8, AR(1))",  dgp_a8, n_reps = 100, nboots = 1000)
    t21 <- run_width_parity(n_reps = 50, nboots = 500)

    out_rds <- file.path(out_dir, sprintf("coverage_para_error_%s.rds", ts))
    saveRDS(list(t19 = t19, t20 = t20, t21 = t21), out_rds)

    ## Long-format CSV summary for diff/audit
    rows <- list()
    for (res in list(t19, t20)) {
        thr_low <- if (grepl("^T19", res$label)) 0.90 else 0.91
        for (pm in rownames(res$cov_mat)) for (m in colnames(res$cov_mat)) {
            rows[[length(rows) + 1L]] <- data.frame(
                test = res$label, dgp = sub(".*\\(([^)]+)\\).*", "\\1", res$label),
                para.error = pm, ci.method = m,
                coverage = res$cov_mat[pm, m], width = res$width_mat[pm, m],
                threshold_low = thr_low, threshold_high = 0.99,
                pass = res$cov_mat[pm, m] >= thr_low & res$cov_mat[pm, m] <= 0.99,
                stringsAsFactors = FALSE)
        }
    }
    summary_df <- do.call(rbind, rows)
    csv_path <- file.path(out_dir, sprintf("coverage_para_error_%s.csv", ts))
    write.csv(summary_df, csv_path, row.names = FALSE)

    cat("\n", strrep("=", 60), "\n", sep = "")
    cat("OVERALL: T19 ", t19$pass, "  T20 ", t20$pass, "  T21 ", t21$pass, "\n", sep = "")
    cat("CSV:  ", csv_path, "\n", sep = "")
    cat("RDS:  ", out_rds, "\n", sep = "")
    cat(strrep("=", 60), "\n", sep = "")
    invisible(list(t19 = t19, t20 = t20, t21 = t21,
                   csv = csv_path, rds = out_rds))
}

## Run when sourced as a script (not when load_all'd interactively)
if (sys.nframe() == 0L || identical(environmentName(parent.frame()), "R_GlobalEnv")) {
    run_all()
}
