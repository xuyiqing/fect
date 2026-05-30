## ============================================================================
## Phase B / Phase B-CFE: partial warm-start validation (v2.4.3)
##
## Implementation: warm.start = c("none", "linear") on fect().
## "linear" enables partial warm-start in the bootstrap loop:
##   - per-replicate fits seed the EM with the auxiliary-only prediction
##     surface from the cached main fit (Y.ct - factor %*% t(lambda));
##   - factor pair (F, Lambda) cold-starts via fresh per-replicate SVD
##     on the auxiliary-corrected residual.
##
## Phase A (full warm-start, anchoring F, Lambda from main fit) FAILED
## empirical variance-neutrality (cold-vs-warm SE diverged 45-250%
## relative). Partial warm-start preserves the basin re-randomization
## (because Y_b varies per replicate) while still warming all the
## auxiliaries that are deterministic functions of (F, Lambda) given Y.
##
## Pass criteria (per design):
##   - point estimate (att.avg, est.att[, "ATT"]) byte-identical
##   - max relative S.E. diff < 5%
##   - speedup (cold_time / warm_time) >= 2x
##
## RUN WHEN: validating partial warm-start. ~30 min total wall time.
## ============================================================================

suppressPackageStartupMessages({ library(devtools) })
setwd("/Users/xyq/GitHub/fect-warmstart")
devtools::load_all(quiet = TRUE)

OUT_DIR <- "/tmp/fect-partial-warm"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------ Utilities --

compare_fits <- function(fit_cold, fit_warm, cold_time, warm_time, label) {
    diff_att_avg  <- abs(fit_cold$att.avg - fit_warm$att.avg)
    diff_est_att  <- max(abs(fit_cold$est.att[, "ATT"] -
                             fit_warm$est.att[, "ATT"]))
    diff_eff_boot <- if (!is.null(fit_cold$eff.boot) &&
                          !is.null(fit_warm$eff.boot))
        max(abs(fit_cold$eff.boot - fit_warm$eff.boot), na.rm = TRUE) else NA_real_
    se_cold <- fit_cold$est.att[, "S.E."]
    se_warm <- fit_warm$est.att[, "S.E."]
    se_rel_diff <- max(abs(se_cold - se_warm) / pmax(se_cold, 1e-10))
    speedup <- cold_time / warm_time

    pass_pt   <- diff_att_avg < 1e-8 && diff_est_att < 1e-8
    pass_se   <- se_rel_diff < 0.05
    pass_speed <- speedup >= 2

    res <- list(
        label          = label,
        diff_att_avg   = diff_att_avg,
        diff_est_att   = diff_est_att,
        diff_eff_boot  = diff_eff_boot,
        se_rel_diff    = se_rel_diff,
        speedup        = speedup,
        cold_time      = cold_time,
        warm_time      = warm_time,
        pass_pt        = pass_pt,
        pass_se        = pass_se,
        pass_speed     = pass_speed,
        pass_overall   = pass_pt && pass_se && pass_speed
    )

    cat(sprintf("\n[%s]\n", label))
    cat(sprintf("  cold time: %.2fs   warm time: %.2fs   speedup: %.2fx %s\n",
                cold_time, warm_time, speedup,
                if (pass_speed) "PASS" else "FAIL"))
    cat(sprintf("  diff att.avg:    %.2e   %s\n",
                diff_att_avg, if (pass_pt) "PASS" else "FAIL"))
    cat(sprintf("  diff est.att:    %.2e\n", diff_est_att))
    cat(sprintf("  diff eff.boot:   %.4e\n", diff_eff_boot))
    cat(sprintf("  rel S.E. diff:   %.2f%% %s\n",
                100 * se_rel_diff, if (pass_se) "PASS" else "FAIL"))
    cat(sprintf("  OVERALL:         %s\n",
                if (res$pass_overall) "PASS" else "FAIL"))
    res
}

## ------------------------------------------------------------ Phase B (IFE) --

run_phase_b_ife <- function(n_reps = 30, seed = 42) {
    cat("\n========== Phase B: IFE r=2 on simdata, n=", n_reps, " ==========\n",
        sep = "")
    data(simdata)

    ## Use tol = 1e-5 for both cold and warm so the bootstrap distributions
    ## estimate the same target (full convergence). At fect's default
    ## tol = 1e-3 the EM stops basin-dependent and warm-vs-cold are
    ## different estimators (see Phase B finding 2026-05-02).
    common_args <- list(
        formula = Y ~ D + X1 + X2, data = simdata, index = c("id","time"),
        method = "ife", r = 2, force = "two-way",
        se = TRUE, vartype = "bootstrap", nboots = n_reps,
        parallel = FALSE, keep.sims = TRUE, CV = FALSE,
        tol = 1e-5, max.iteration = 5000
    )

    set.seed(seed)
    t0 <- Sys.time()
    fit_cold <- do.call(fect, common_args)
    cold_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    set.seed(seed)
    t0 <- Sys.time()
    fit_warm <- do.call(fect, c(common_args, list(warm.start = "linear")))
    warm_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    compare_fits(fit_cold, fit_warm, cold_time, warm_time,
                 label = "Phase B IFE simdata r=2")
}

## ------------------------------------------------------------ Phase B-CFE --

run_phase_b_cfe <- function(dgp_name, n_reps = 30, seed = 42) {
    cat("\n========== Phase B-CFE: ", dgp_name, ", n=", n_reps, " ==========\n",
        sep = "")

    dgp_specs <- list(
        sim_region = list(
            data_name = "sim_region",
            formula   = Y ~ D,
            sfe       = "region",
            cfe       = NULL,
            Q         = NULL,
            r         = 0,
            note      = "Multi-level FE (region); no factors"
        ),
        sim_trend = list(
            data_name = "sim_trend",
            formula   = Y ~ D,
            sfe       = NULL,
            cfe       = list(c("id","year")),
            Q         = NULL,
            r         = 0,
            note      = "Q-heavy stand-in: cfe=list(c(id,year)) demonstrates kappa path"
        ),
        sim_linear = list(
            data_name = "sim_linear",
            formula   = Y ~ D,
            sfe       = NULL,
            cfe       = list(c("id", "year")),
            Q         = NULL,
            r         = 0,
            note      = "Q-light stand-in: linear-trend approximated via cfe"
        ),
        simdata = list(
            data_name = "simdata",
            formula   = Y ~ D + X1 + X2,
            sfe       = NULL,
            cfe       = NULL,
            Q         = NULL,
            r         = 2,
            note      = "Factor-only DGP with covariates"
        ),
        sim_gsynth = list(
            data_name = "sim_gsynth",
            formula   = Y ~ D + X1 + X2,
            sfe       = NULL,
            cfe       = NULL,
            Q         = NULL,
            r         = 2,
            method_override = "gsynth",
            note      = "Gsynth path (fect_nevertreated)"
        )
    )

    spec <- dgp_specs[[dgp_name]]
    if (is.null(spec)) stop("Unknown DGP: ", dgp_name)
    cat("  spec: ", spec$note, "\n", sep = "")

    e <- new.env()
    do.call(data, list(spec$data_name, package = "fect", envir = e))
    df <- get(spec$data_name, envir = e)

    method_use <- if (!is.null(spec$method_override)) spec$method_override else "cfe"

    fit_args <- list(
        formula  = spec$formula,
        data     = df,
        index    = c("id", "time"),
        method   = method_use,
        force    = "two-way",
        sfe      = spec$sfe,
        cfe      = spec$cfe,
        Q        = spec$Q,
        r        = spec$r,
        se       = TRUE,
        vartype  = "bootstrap",
        nboots   = n_reps,
        parallel = FALSE,
        keep.sims = TRUE,
        CV       = FALSE,
        tol      = 1e-5,
        max.iteration = 5000
    )

    cold_args <- c(fit_args, list(warm.start = "none"))
    warm_args <- c(fit_args, list(warm.start = "linear"))

    set.seed(seed)
    t0 <- Sys.time()
    fit_cold <- tryCatch(do.call(fect, cold_args), error = function(e) e)
    if (inherits(fit_cold, "error")) {
        cat("  cold-start ERROR:", conditionMessage(fit_cold), "\n")
        return(list(label = paste0("Phase B-CFE ", dgp_name),
                    error = conditionMessage(fit_cold), pass_overall = FALSE))
    }
    cold_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    set.seed(seed)
    t0 <- Sys.time()
    fit_warm <- tryCatch(do.call(fect, warm_args), error = function(e) e)
    if (inherits(fit_warm, "error")) {
        cat("  warm-start ERROR:", conditionMessage(fit_warm), "\n")
        return(list(label = paste0("Phase B-CFE ", dgp_name),
                    error = conditionMessage(fit_warm), pass_overall = FALSE))
    }
    warm_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    compare_fits(fit_cold, fit_warm, cold_time, warm_time,
                 label = paste0("Phase B-CFE ", dgp_name))
}

## ------------------------------------------------------------ Driver --

run_all <- function(out_dir = OUT_DIR, n_reps = 30) {
    ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
    cat("\n=== Partial warm-start validation (v2.4.3 candidate) ===\n")
    cat("Started: ", as.character(Sys.time()), "\n")
    cat("nboots per cell: ", n_reps, "\n")

    results <- list(
        phase_b_ife = run_phase_b_ife(n_reps = n_reps),
        phase_b_cfe_sim_region = run_phase_b_cfe("sim_region", n_reps = n_reps),
        phase_b_cfe_sim_trend  = run_phase_b_cfe("sim_trend",  n_reps = n_reps),
        phase_b_cfe_sim_linear = run_phase_b_cfe("sim_linear", n_reps = n_reps),
        phase_b_cfe_simdata    = run_phase_b_cfe("simdata",    n_reps = n_reps),
        phase_b_cfe_sim_gsynth = run_phase_b_cfe("sim_gsynth", n_reps = n_reps)
    )

    cat("\n========================================\n")
    cat("Summary\n")
    cat("========================================\n")
    cat(sprintf("%-32s %-7s %-9s %-7s %-9s %-7s %-7s\n",
                "label", "att.diff", "se.rel", "speed", "pass.pt", "pass.se", "pass.spd"))
    for (r in results) {
        if (!is.null(r$error)) {
            cat(sprintf("%-32s ERROR: %s\n", r$label, r$error))
        } else {
            cat(sprintf("%-32s %.2e  %5.2f%%   %4.2fx %5s     %5s     %5s\n",
                        r$label, r$diff_att_avg, 100 * r$se_rel_diff,
                        r$speedup,
                        if (r$pass_pt) "PASS" else "FAIL",
                        if (r$pass_se) "PASS" else "FAIL",
                        if (r$pass_speed) "PASS" else "FAIL"))
        }
    }

    out_file <- file.path(out_dir, sprintf("phase_b_results_%s.rds", ts))
    saveRDS(results, out_file)
    cat("Saved: ", out_file, "\n", sep = "")
    invisible(results)
}

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    n_reps <- if (length(args) > 0) as.integer(args[[1]]) else 100
    run_all(n_reps = n_reps)
}
