## ============================================================================
## Conditional follow-up to run_minimal_coverage.R: re-run any scenario whose
## tail-CI methods (basic / percentile / bc / bca) come in below threshold
## (default 0.93) at nboots = 1000 instead of 200.
##
## This complements run_minimal_coverage.R per the v2.4.2 .check_tail_ci_replicates
## warning: tail-quantile CIs at B = 200 can be erratic (E&T 1987 §3, DiCiccio &
## Efron 1996 §4 recommend B >= 1000).  Normal CI is unaffected by B; jackknife
## is normal-only by construction.  So the rerun targets only A/B/C1 cells with
## tail methods that came in low.
##
## Reads the most recent minimal_summary_K*_nb*_*.csv and dispatches reruns.
## ============================================================================

.script_path_self <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0L) {
        return(normalizePath(sub("^--file=", "", file_arg[1])))
    }
    of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(of)) return(normalizePath(of))
    NA_character_
}
.this <- .script_path_self()
.this_dir <- if (!is.na(.this)) dirname(.this) else getwd()
options(.fect_minimal_no_autorun = TRUE)   # prevent run_all() inside source
source(file.path(.this_dir, "run_minimal_coverage.R"), local = FALSE)
options(.fect_minimal_no_autorun = NULL)
## run_all() at the bottom of run_minimal_coverage.R is gated by sys.nframe()
## so it won't trigger when sourced from inside this script.

read_summary <- function(out_dir = "/tmp/fect-coverage-study") {
    files <- list.files(out_dir, pattern = "^minimal_summary_K.*nb200.*\\.csv$",
                         full.names = TRUE)
    if (length(files) == 0L) {
        stop("No minimal_summary_*_nb200_*.csv found; run run_minimal_coverage.R first.")
    }
    latest <- files[order(file.mtime(files), decreasing = TRUE)][1]
    cat("Reading summary:", latest, "\n")
    read.csv(latest, stringsAsFactors = FALSE)
}

decide_reruns <- function(summary, threshold = 0.93,
                           tail_methods = c("basic", "percentile", "bc", "bca")) {
    fail <- summary$ci.method %in% tail_methods & summary$coverage < threshold
    if (!any(fail)) return(character(0L))
    fail_rows <- summary[fail, ]
    cat("\nFailing tail-CI cells (coverage < ", threshold, "):\n", sep = "")
    print(fail_rows[, c("scenario", "ci.method", "coverage", "mc_se")],
          row.names = FALSE)
    unique(fail_rows$scenario)
}

main <- function(K = 200, nboots_tail = 1000, workers = 16,
                  threshold = 0.93,
                  out_dir = "/tmp/fect-coverage-study") {
    summary <- read_summary(out_dir)
    cat("\n=== Initial summary (nboots = 200) ===\n")
    print(summary[, c("scenario", "ci.method", "coverage", "mc_se",
                       "mean_width")], row.names = FALSE, digits = 4)

    rerun_scenarios <- decide_reruns(summary, threshold = threshold)
    if (length(rerun_scenarios) == 0L) {
        cat("\nAll tail-CI cells >= ", threshold, ". No reruns needed.\n", sep = "")
        return(invisible(summary))
    }
    cat("\nScenarios requiring rerun at nboots = ", nboots_tail, ": ",
        paste(rerun_scenarios, collapse = ", "), "\n", sep = "")

    future::plan(future::multisession, workers = workers)
    options(future.globals.maxSize = 2e9)
    on.exit(future::plan(future::sequential), add = TRUE)

    ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
    rerun_summaries <- list()
    for (scen in rerun_scenarios) {
        seed_base <- switch(scen,
                             "A"      = 1000L,
                             "B"      = 2000L,
                             "C_boot" = 3000L,
                             "C_jack" = 4000L,
                             stop("Unknown scenario: ", scen))
        ## Use the SAME seeds as the original run -- ensures the only changing
        ## factor is nboots (paired comparison).
        label <- switch(scen,
                          "A"      = "A: factor IID parametric",
                          "B"      = "B: factor AR(1) rho=0.8 parametric",
                          "C_boot" = "C1: TWFE AR(1) rho=0.5 bootstrap",
                          "C_jack" = "C2: TWFE AR(1) rho=0.5 jackknife",
                          paste(scen, "(rerun)"))
        result <- run_scenario(scen, paste0(label, "  [nboots=", nboots_tail, "]"),
                                K = K, nboots = nboots_tail, workers = workers,
                                base_seed = seed_base,
                                ci_methods = if (scen == "C_jack") "normal"
                                             else CI_METHODS_FULL)
        df_csv <- file.path(out_dir,
                             sprintf("minimal_%s_K%d_nb%d_%s.csv",
                                     scen, K, nboots_tail, ts))
        write.csv(result$df, df_csv, row.names = FALSE)
        cat("  wrote", df_csv, "\n")
        rerun_summaries[[scen]] <- result$summary
    }

    rerun_df <- do.call(rbind, rerun_summaries)
    sum_csv <- file.path(out_dir,
                         sprintf("minimal_tail_rerun_summary_K%d_nb%d_%s.csv",
                                 K, nboots_tail, ts))
    write.csv(rerun_df, sum_csv, row.names = FALSE)

    cat("\n", strrep("=", 78), "\n", sep = "")
    cat("RERUN SUMMARY  (nboots = ", nboots_tail, ")\n", sep = "")
    cat(strrep("=", 78), "\n", sep = "")
    print(rerun_df[, c("scenario", "ci.method", "coverage", "mc_se",
                        "mean_width")], row.names = FALSE, digits = 4)
    cat(sprintf("\nRerun summary CSV: %s\n", sum_csv))
    cat(strrep("=", 78), "\n", sep = "")

    invisible(rerun_df)
}

if (sys.nframe() == 0L || identical(environmentName(parent.frame()), "R_GlobalEnv")) {
    main()
}
