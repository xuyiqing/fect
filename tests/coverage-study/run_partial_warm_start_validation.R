## ============================================================================
## Phase B / Phase B-CFE: partial warm-start validation
##
## Tests the design where auxiliary parameters are warm-started from the main
## fit (alpha, xi, beta, [gamma, kappa, multi-FE for CFE], missing-cell fillers)
## but the factor pair (F, Lambda) cold-starts via fresh SVD per replicate.
##
## Premise (post-Phase-A): full warm-start (carrying main fit's F, Lambda)
## broke variance-neutrality of the bootstrap because all replicates landed
## in the same factor basin, deflating SE estimates by 45-250% relative.
## Partial warm-start preserves per-replicate basin variability while still
## getting most of the speedup, by NOT carrying (F, Lambda).
##
## Status: Phase B research script. Implementation of warm.start in fect() is
## NOT yet shipped; this script prototypes the idea via direct manipulation
## of the EM input. If Phase B passes, ship in v2.4.3.
##
## Design doc: statsclaw-workspace/fect/ref/v242-warm-start-investigation/
##             partial-warm-design.md
##
## RUN WHEN: investigating whether partial warm-start preserves variance.
## DO NOT RUN ON: routine devtools::test().
## Wall time: Phase B IFE ~5 min; Phase B-CFE ~30 min across 5 DGPs.
##
## Output: /tmp/fect-partial-warm/results_<timestamp>.rds + stdout summary.
## ============================================================================

suppressPackageStartupMessages({ library(devtools) })
setwd("/Users/xyq/GitHub/fect")
devtools::load_all(quiet = TRUE)

OUT_DIR <- "/tmp/fect-partial-warm"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------ Utilities --

## Construct the auxiliary-only prediction surface from a fitted model.
## This is what we pass as fit_init for partial warm-start: it contains
## alpha + xi + (X*beta + Z*gamma + Q*kappa + multi_FE if applicable) +
## main-fit's Y0_hat at missing cells, but NOT the F*Lambda factor product.
##
## When the EM seeds with this, the first SVD step finds fresh per-replicate
## (F, Lambda) on the residual Y_b - aux_surface.
build_aux_surface <- function(fit) {
    ## Reconstruct alpha + xi from fit
    ## fit$alpha.tr / fit$xi: vectors per unit / time
    ## fit$beta:             covariate coefficients
    ## fit$Y.ct (or similar): full counterfactual surface = aux + F*Lambda
    ## fit$factor / fit$lambda.tr: factor matrices
    ##
    ## Aux surface = Y.ct - F %*% t(Lambda)   (subtract factor part)

    if (is.null(fit$Y.ct)) {
        stop("fit does not expose Y.ct; cannot build aux surface")
    }
    aux <- fit$Y.ct
    if (!is.null(fit$factor) && !is.null(fit$lambda.tr) &&
        ncol(fit$factor) > 0) {
        ## Factor product = F %*% t(Lambda)
        aux <- aux - fit$factor %*% t(fit$lambda.tr)
    }
    aux
}

## Apply partial warm-start to a fect() call by passing fit_init = aux surface.
## Currently fect() does not expose fit_init in its public API; this prototype
## requires a small internal patch or a direct call to the C++ entry. For now
## we ASSUME fect() has a hidden `.fit_init` argument that gets routed to
## inter_fe_ub. The actual production rollout will expose this as
## `warm.start = "linear"` in fect().
fect_partial_warm <- function(formula, data, index, aux_init, ...) {
    ## Placeholder: in production this would call:
    ## fect(..., warm.start = "linear", .aux_init = aux_init)
    ##
    ## For the prototype, we simulate partial warm-start by externally
    ## subtracting the auxiliary surface from Y before fitting. This is
    ## NOT identical to true partial warm-start (because the EM won't
    ## refine the auxiliaries internally), but it captures the key
    ## statistical question: does (F, Lambda) basin vary per replicate
    ## when initialized from the auxiliary-corrected matrix?

    stop(paste(
        "fect_partial_warm requires either",
        "(a) a fect() patch exposing fit_init via .aux_init, or",
        "(b) direct call to inter_fe_ub with the auxiliary surface.",
        "Implement before running Phase B; see partial-warm-design.md"
    ))
}

## ------------------------------------------------------------ Phase B (IFE) --

run_phase_b_ife <- function(n_reps = 30, seed = 42) {
    cat("\n=== Phase B: IFE r=2 on simdata, ", n_reps, " reps ===\n", sep = "")

    data(simdata, package = "fect")

    ## Main fit (cold-start, baseline)
    cat("Main fit (cold)...\n")
    set.seed(seed)
    t0 <- Sys.time()
    fit_main <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                     method = "ife", r = 2, force = "two-way",
                     se = TRUE, vartype = "bootstrap",
                     nboots = n_reps, parallel = FALSE,
                     keep.sims = TRUE, CV = FALSE)
    main_cold_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  cold-start wall: %.1f sec\n", main_cold_time))

    ## Build auxiliary surface for partial warm-start
    aux_surface <- build_aux_surface(fit_main)
    cat(sprintf("  aux_surface dim: %s, range [%.2f, %.2f]\n",
                paste(dim(aux_surface), collapse = " x "),
                min(aux_surface, na.rm = TRUE), max(aux_surface, na.rm = TRUE)))

    ## Partial-warm fit (REQUIRES IMPLEMENTATION)
    cat("Partial-warm fit (NOT YET IMPLEMENTED)...\n")
    cat("  See fect_partial_warm() docstring; requires fect() patch.\n")
    ## set.seed(seed)
    ## t0 <- Sys.time()
    ## fit_warm <- fect_partial_warm(Y ~ D + X1 + X2, data = simdata,
    ##                                index = c("id", "time"),
    ##                                aux_init = aux_surface,
    ##                                method = "ife", r = 2, force = "two-way",
    ##                                se = TRUE, vartype = "bootstrap",
    ##                                nboots = n_reps, parallel = FALSE,
    ##                                keep.sims = TRUE, CV = FALSE)
    ## warm_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    ## Comparison metrics (when fit_warm is available)
    ## diff_att_avg <- abs(fit_main$att.avg - fit_warm$att.avg)
    ## diff_est_att <- max(abs(fit_main$est.att[, "ATT"] - fit_warm$est.att[, "ATT"]))
    ## diff_eff_boot <- max(abs(fit_main$eff.boot - fit_warm$eff.boot))
    ## se_relative_diff <- max(abs(fit_main$est.att[, "S.E."] - fit_warm$est.att[, "S.E."]) /
    ##                         pmax(fit_main$est.att[, "S.E."], 1e-10))
    ## speedup <- main_cold_time / warm_time

    list(main_cold_time = main_cold_time,
         aux_surface = aux_surface,
         status = "implementation_pending")
}

## ------------------------------------------------------------ Phase B-CFE --

run_phase_b_cfe <- function(dgp_name, n_reps = 30, seed = 42) {
    cat("\n=== Phase B-CFE: ", dgp_name, ", ", n_reps, " reps ===\n", sep = "")

    ## DGP -> CFE specification mapping
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
            cfe       = NULL,
            Q         = "id",     ## sinusoidal trend basis on id
            r         = 0,
            note      = "Q-heavy: unit-specific sinusoidal trends"
        ),
        sim_linear = list(
            data_name = "sim_linear",
            formula   = Y ~ D,
            sfe       = NULL,
            cfe       = NULL,
            Q         = "id",
            Q.type    = "linear",
            r         = 0,
            note      = "Q-light: unit-specific linear trends"
        ),
        simdata = list(
            data_name = "simdata",
            formula   = Y ~ D + X1 + X2,
            sfe       = NULL,
            cfe       = NULL,
            Q         = NULL,
            r         = 2,
            note      = "Factor-only DGP; CFE without Q (equivalent to IFE)"
        ),
        sim_gsynth = list(
            data_name = "sim_gsynth",
            formula   = Y ~ D + X1 + X2,
            sfe       = NULL,
            cfe       = NULL,
            Q         = NULL,
            r         = 2,
            time.component.from = "nevertreated",
            note      = "Factor + covariates; gsynth-style"
        )
    )

    spec <- dgp_specs[[dgp_name]]
    if (is.null(spec)) stop("Unknown DGP: ", dgp_name)

    cat("  spec: ", spec$note, "\n", sep = "")

    e <- new.env()
    do.call(data, list(spec$data_name, package = "fect", envir = e))
    df <- get(spec$data_name, envir = e)

    fit_args <- list(
        formula  = spec$formula,
        data     = df,
        index    = c("id", "time"),
        method   = "cfe",
        force    = "two-way",
        sfe      = spec$sfe,
        cfe      = spec$cfe,
        Q        = spec$Q,
        Q.type   = spec$Q.type,
        r        = spec$r,
        se       = TRUE,
        vartype  = "bootstrap",
        nboots   = n_reps,
        parallel = FALSE,
        keep.sims = TRUE,
        CV       = FALSE
    )
    if (!is.null(spec$time.component.from)) {
        fit_args$time.component.from <- spec$time.component.from
    }

    cat("Cold-start fit...\n")
    set.seed(seed)
    t0 <- Sys.time()
    fit_cold <- do.call(fect, fit_args)
    cold_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  cold-start wall: %.1f sec\n", cold_time))

    cat("Partial-warm fit (NOT YET IMPLEMENTED)...\n")
    ## See run_phase_b_ife for the same comment.

    list(dgp = dgp_name, cold_time = cold_time, status = "implementation_pending")
}

## ------------------------------------------------------------ Driver --

run_all <- function(out_dir = OUT_DIR) {
    ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
    cat("Phase B / Phase B-CFE prototype run starting at ", as.character(Sys.time()), "\n", sep = "")
    cat("This script is a SCAFFOLD --- it requires fect() to expose fit_init\n")
    cat("via a `warm.start` argument or `.aux_init` parameter before the\n")
    cat("partial-warm comparison can run. See:\n")
    cat("  statsclaw-workspace/fect/ref/v242-warm-start-investigation/\n")
    cat("    partial-warm-design.md\n\n")

    results <- list(
        phase_b_ife = run_phase_b_ife(),
        phase_b_cfe_sim_region = run_phase_b_cfe("sim_region"),
        phase_b_cfe_sim_trend  = run_phase_b_cfe("sim_trend"),
        phase_b_cfe_sim_linear = run_phase_b_cfe("sim_linear"),
        phase_b_cfe_simdata    = run_phase_b_cfe("simdata"),
        phase_b_cfe_sim_gsynth = run_phase_b_cfe("sim_gsynth")
    )

    out_file <- file.path(out_dir, sprintf("phase_b_results_%s.rds", ts))
    saveRDS(results, out_file)
    cat("Saved: ", out_file, "\n", sep = "")

    invisible(results)
}

## Run when sourced as a script
if (sys.nframe() == 0L) {
    run_all()
}
