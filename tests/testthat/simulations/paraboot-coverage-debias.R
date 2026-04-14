###############################################################################
## Parametric Bootstrap Debiasing Coverage Simulation — Stage 4 Harness
## Run: 2026-04-14-stage4-debiasing-poc
## Workflow: 11 (code + simulation — builder || simulator, then tester)
##
## Purpose: Empirical 95% CI coverage under parametric bootstrap with and
##          without sample-splitting debiasing (split_residuals=TRUE/FALSE).
##          8-cell focused grid (D1-D8).
##
## Estimator interface: fect(..., split_residuals = TRUE/FALSE)
##   split_residuals = FALSE  => standard parametric bootstrap (baseline)
##   split_residuals = TRUE   => debiased: control units split into two halves;
##                              half-A for factor estimation, half-B for
##                              residual sampling (builder provides this flag)
##
## Usage (from repo root):
##   Rscript tests/testthat/simulations/paraboot-coverage-debias.R \
##     [--fect-path <path>] [--cores <n>] [--out-dir <dir>] \
##     [--n-sim <n>] [--only-cells <D1,D2,...>]
##
## Arguments (all optional):
##   --fect-path  : path to fect package root (default: ".")
##   --cores      : outer mclapply workers (default: 14)
##   --out-dir    : directory for RDS outputs (default: same dir as this script)
##   --n-sim      : override N_sim for all cells (default: per grid spec)
##   --only-cells : comma-separated cell IDs to run, e.g. "D2,D4,D6"
##
## Per-cell checkpointing:
##   Each completed cell is saved atomically as <out-dir>/cell-<cell_id>.rds.
##   On re-run, existing checkpoints are loaded and skipped (no re-run).
##   Use --only-cells to target specific cells in a restart.
##
## Gate C handling (D1, D3, D5):
##   Cells D1, D3, D5 have split_residuals=FALSE and tcf="notyettreated".
##   With post-Stage-4 code, Gate C blocks these calls. The harness wraps
##   each fect() call in tryCatch and records gate_blocked=TRUE if a Gate C
##   error is encountered. Stage 3 reference coverage values are hard-coded.
##
## Reproducibility:
##   Master seed: 20260414. Cell seed = master + (cell_index-1)*10000.
##   Rep seed = cell_seed + rep_index.  RNG: Mersenne-Twister (R default).
##
## Output files:
##   <out-dir>/cell-D<n>.rds            — per-cell checkpoint (raw + summary)
##   <out-dir>/debias-results.rds       — aggregated raw results
##   <out-dir>/debias-summary.rds       — summary data.frame (8 rows)
##
## Do NOT run this script directly — tester invokes it after builder+simulator
## merge back. See test-spec.md for the full validation protocol.
###############################################################################

## ---- Smoke test helper -------------------------------------------------------
## smoke_test() is called by tester to verify the API contract before the full
## run. Runs R=5 replications on D2 (ife+nt, split=TRUE, N=50) using a
## temporary output directory. Returns invisibly; prints a summary line.
##
## Invoke:
##   source("tests/testthat/simulations/paraboot-coverage-debias.R")
##   smoke_test(fect_path = ".", cores = 1L)
##
## Expected: no errors, valid numeric coverage in [0, 1], non-NULL att_est.
smoke_test <- function(fect_path = ".", cores = 1L) {
  message("=== SMOKE TEST: D2 (ife+nt+split=TRUE, N=50, R=5) ===")
  td <- tempfile(pattern = "debias_smoke_")
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  args_save <- commandArgs
  ## We run the cell directly rather than re-sourcing, to avoid re-parsing args.
  ## Load dependencies first.
  if (!exists("generate_ife_dgp", mode = "function")) {
    dgp_file_smoke <- file.path(fect_path,
                                "tests/testthat/simulations/dgp-helpers.R")
    if (!file.exists(dgp_file_smoke)) {
      dgp_file_smoke <- "tests/testthat/simulations/dgp-helpers.R"
    }
    source(dgp_file_smoke)
  }
  if (!exists("fect", mode = "function")) {
    suppressMessages(devtools::load_all(fect_path))
  }

  SMOKE_SEED <- 20270414L  # D2 cell seed
  SMOKE_TAU  <- 1.0

  reps_smoke <- parallel::mclapply(seq_len(5L), function(s) {
    set.seed(SMOKE_SEED + s)
    simdf <- generate_ife_dgp(N = 50L, TT = 20L, Ntr = 10L,
                               tau = SMOKE_TAU, sigma_eps = 1.0,
                               r = 2L, seed = NULL)
    out <- tryCatch(
      suppressWarnings(suppressMessages(fect(
        Y ~ D,
        data    = simdf,
        index   = c("id", "time"),
        method  = "ife",
        r       = 2L,
        CV      = FALSE,
        force   = "two-way",
        se      = TRUE,
        vartype = "parametric",
        nboots  = 10L,          # tiny B for smoke test
        parallel = FALSE,
        split_residuals = TRUE
      ))),
      error = function(e)
        structure(conditionMessage(e), class = "try-error", condition = e)
    )
    extract_result_debias(out, tau_true = SMOKE_TAU)
  }, mc.cores = cores, mc.set.seed = FALSE)

  n_fail  <- sum(sapply(reps_smoke, function(r) isTRUE(r$failed)))
  n_gate  <- sum(sapply(reps_smoke, function(r) isTRUE(r$gate_blocked)))
  att_ok  <- !all(is.na(sapply(reps_smoke, function(r) r$att_est)))

  message(sprintf("  reps=5  failures=%d  gate_blocked=%d  att_non-NA=%s",
                  n_fail, n_gate, if (att_ok) "YES" else "NO"))

  if (n_fail > 0) {
    msgs <- unique(sapply(reps_smoke, function(r) r$fail_msg))
    msgs <- msgs[!is.na(msgs)]
    message("  Failure messages: ", paste(msgs[seq_len(min(3L, length(msgs)))],
                                          collapse = " | "))
  }

  invisible(reps_smoke)
}

## ---- Parse command-line arguments -------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(args, flag, default) {
  idx <- which(args == flag)
  if (length(idx) == 1L && idx < length(args)) args[idx + 1L] else default
}

fect_path       <- get_arg(args, "--fect-path", ".")
n_cores         <- as.integer(get_arg(args, "--cores",   "14"))
out_dir         <- get_arg(args, "--out-dir",  ".")
n_sim_override  <- get_arg(args, "--n-sim",    NA_character_)
only_cells_str  <- get_arg(args, "--only-cells", NA_character_)

# Parse --only-cells into a character vector (NULL = run all)
only_cells <- if (!is.na(only_cells_str)) {
  trimws(strsplit(only_cells_str, ",")[[1]])
} else {
  NULL
}

if (!is.null(only_cells)) {
  message(sprintf("--only-cells filter active: %s",
                  paste(only_cells, collapse = ", ")))
}

## ---- Ensure output directory exists -----------------------------------------
dir.create(out_dir,                     recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tmp"),   recursive = TRUE, showWarnings = FALSE)

## ---- Load package from source ------------------------------------------------
message("Loading fect from: ", normalizePath(fect_path))
suppressMessages(devtools::load_all(fect_path))

# Emit the fect git commit for the environment record
fect_git_commit <- tryCatch(
  system2("git", c("-C", fect_path, "log", "-1", "--oneline"),
          stdout = TRUE, stderr = FALSE),
  error = function(e) "UNKNOWN"
)
message("fect git HEAD: ", paste(fect_git_commit, collapse = " "))

## ---- Source DGP helpers ------------------------------------------------------
dgp_file <- file.path(fect_path, "tests/testthat/simulations/dgp-helpers.R")
if (!file.exists(dgp_file)) {
  dgp_file <- "tests/testthat/simulations/dgp-helpers.R"
}
if (!file.exists(dgp_file)) {
  stop("dgp-helpers.R not found. Looked in: ", dgp_file)
}
message("Sourcing DGP helpers from: ", dgp_file)
source(dgp_file)

## ---- Master seed & seed derivation -------------------------------------------
MASTER_SEED <- 20260414L

# Cell index: D1=1, D2=2, ..., D8=8
derive_cell_seed <- function(cell_index) {
  MASTER_SEED + (cell_index - 1L) * 10000L
}

derive_rep_seed <- function(cell_seed, rep_index) {
  cell_seed + rep_index
}

## ---- Stage 3 reference constants (for Gate C blocked cells) ------------------
## D1 (ife+nt, N=50, split=FALSE): from Stage 3 simulation.md row 9
## D3 (ife+nt, N=100, split=FALSE): Stage 3 did not include this cell
## D5 (cfe+nt, N=50, split=FALSE): Stage 3 gated before running; no result
##
## Source: runs/2026-04-13-paraboot-coverage/simulation.md
STAGE3_REF <- list(
  D1 = list(
    coverage_rate = 0.806,
    N_sim         = 500L,
    tau           = 1.0,
    note          = "Stage3 C3_ife_nt_N50_tau1.0 — row 9 of simulation.md"
  ),
  D3 = list(
    coverage_rate = NA_real_,
    N_sim         = 0L,
    tau           = 1.0,
    note          = "Stage3 did not run ife+nt at N=100 (Gate C decision made on N=50 evidence)"
  ),
  D5 = list(
    coverage_rate = NA_real_,
    N_sim         = 0L,
    tau           = 1.0,
    note          = "Stage3 gated cfe+notyet before running; no empirical baseline available"
  )
)

## ---- Checkpoint helpers -------------------------------------------------------
cell_rds_path <- function(cell_id, out_dir) {
  file.path(out_dir, paste0("cell-", cell_id, ".rds"))
}

# Atomic write: write to tmp file, then rename (avoids partial reads on crash)
save_cell_checkpoint <- function(result_obj, cell_id, out_dir) {
  tmp_dir   <- file.path(out_dir, "tmp")
  tmp_file  <- file.path(tmp_dir, paste0("cell-", cell_id, "-tmp.rds"))
  final_file <- cell_rds_path(cell_id, out_dir)
  saveRDS(result_obj, tmp_file)
  file.rename(tmp_file, final_file)
  invisible(final_file)
}

## ---- extract_result_debias() -------------------------------------------------
## Extracts coverage, ATT, CI, and SE from a fect() output object.
## Returns a named list with consistent fields for both success and failure.
##
## Field names verified against fect() output:
##   out$att.avg          — scalar ATT estimate
##   out$est.avg          — matrix with row 1 = att.avg stats
##     columns: "ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value"
##   out$est.att          — matrix with row "att.avg"
##     columns: "CI.lower", "CI.upper" (spec Section 1 alternative accessor)
##
## NOTE: extract_coverage_att() from dgp-helpers.R uses out$est.avg (Stage 3
## convention). This function adds a fallback to out$est.att["att.avg", ...]
## as specified in sim-spec.md Section 1, and adds gate_blocked tracking.
extract_result_debias <- function(out, tau_true) {
  if (inherits(out, "try-error")) {
    err_msg      <- as.character(out)
    gate_blocked <- grepl("not valid when", err_msg, fixed = FALSE) ||
                    grepl("Gate C",         err_msg, fixed = FALSE)
    return(list(
      covered      = NA_integer_,
      att_est      = NA_real_,
      ci_lower     = NA_real_,
      ci_upper     = NA_real_,
      se_est       = NA_real_,
      failed       = TRUE,
      gate_blocked = gate_blocked,
      fail_msg     = err_msg
    ))
  }

  # --- Primary extraction path: out$est.avg (row 1) ---
  att_est  <- tryCatch(as.numeric(out$att.avg), error = function(e) NA_real_)
  ci_lower <- NA_real_
  ci_upper <- NA_real_
  se_est   <- NA_real_

  if (!is.null(out$est.avg)) {
    avg_mat <- out$est.avg
    if (is.matrix(avg_mat) && nrow(avg_mat) >= 1L) {
      if ("CI.lower" %in% colnames(avg_mat))
        ci_lower <- as.numeric(avg_mat[1L, "CI.lower"])
      if ("CI.upper" %in% colnames(avg_mat))
        ci_upper <- as.numeric(avg_mat[1L, "CI.upper"])
      if ("S.E." %in% colnames(avg_mat))
        se_est   <- as.numeric(avg_mat[1L, "S.E."])
      if (is.na(att_est) && "ATT.avg" %in% colnames(avg_mat))
        att_est  <- as.numeric(avg_mat[1L, "ATT.avg"])
    }
  }

  # --- Fallback: out$est.att["att.avg", ...] (sim-spec Section 1) ---
  if ((is.na(ci_lower) || is.na(ci_upper)) && !is.null(out$est.att)) {
    att_mat <- out$est.att
    if (is.matrix(att_mat) && "att.avg" %in% rownames(att_mat)) {
      if (is.na(ci_lower) && "CI.lower" %in% colnames(att_mat))
        ci_lower <- as.numeric(att_mat["att.avg", "CI.lower"])
      if (is.na(ci_upper) && "CI.upper" %in% colnames(att_mat))
        ci_upper <- as.numeric(att_mat["att.avg", "CI.upper"])
    }
  }

  # --- Derive SE from CI width if S.E. not available ---
  if (is.na(se_est) && !is.na(ci_lower) && !is.na(ci_upper)) {
    se_est <- (ci_upper - ci_lower) / (2.0 * 1.96)
  }

  failed <- is.na(att_est)

  covered <- if (!failed && !is.na(ci_lower) && !is.na(ci_upper)) {
    as.integer(ci_lower <= tau_true & tau_true <= ci_upper)
  } else {
    NA_integer_
  }

  list(
    covered      = covered,
    att_est      = att_est,
    ci_lower     = ci_lower,
    ci_upper     = ci_upper,
    se_est       = se_est,
    failed       = failed || is.na(covered),
    gate_blocked = FALSE,
    fail_msg     = NA_character_
  )
}

## ---- Cell aggregation helper -------------------------------------------------
aggregate_cell_debias <- function(reps_list, tau_true) {
  n_total      <- length(reps_list)
  n_failed     <- sum(sapply(reps_list, function(r) isTRUE(r$failed)),
                      na.rm = TRUE)
  n_gate       <- sum(sapply(reps_list, function(r) isTRUE(r$gate_blocked)),
                      na.rm = TRUE)
  n_valid      <- n_total - n_failed

  att_vals  <- sapply(reps_list, function(r) r$att_est)
  cov_vals  <- sapply(reps_list, function(r) r$covered)
  se_vals   <- sapply(reps_list, function(r) r$se_est)

  att_valid <- att_vals[!is.na(att_vals)]
  cov_valid <- cov_vals[!is.na(cov_vals)]
  se_valid  <- se_vals[!is.na(se_vals)]

  mean_att  <- if (length(att_valid) > 0L) mean(att_valid)    else NA_real_
  sd_att    <- if (length(att_valid) > 1L) sd(att_valid)      else NA_real_
  mean_se   <- if (length(se_valid)  > 0L) mean(se_valid)     else NA_real_
  cov_rate  <- if (length(cov_valid) > 0L) mean(cov_valid)    else NA_real_

  mc_se_cov <- if (!is.na(cov_rate) && length(cov_valid) > 0L) {
    sqrt(cov_rate * (1.0 - cov_rate) / length(cov_valid))
  } else {
    NA_real_
  }

  bias      <- if (!is.na(mean_att)) mean_att - tau_true else NA_real_
  se_ratio  <- if (!is.na(mean_se) && !is.na(sd_att) && sd_att > 0.0) {
    mean_se / sd_att
  } else {
    NA_real_
  }
  fail_rate <- n_failed / n_total
  gate_rate <- n_gate   / n_total

  list(
    n_total   = n_total,
    n_valid   = n_valid,
    n_failed  = n_failed,
    n_gate    = n_gate,
    cov_rate  = cov_rate,
    mc_se_cov = mc_se_cov,
    mean_att  = mean_att,
    bias      = bias,
    sd_att    = sd_att,
    mean_se   = mean_se,
    se_ratio  = se_ratio,
    fail_rate = fail_rate,
    gate_rate = gate_rate
  )
}

## ---- 8-cell scenario grid (D1-D8) -------------------------------------------
## Columns: cell_id, cell_index, method, tcf, dgp_type, N, TT, Ntr, tau,
##          sigma_eps, split_residuals, N_sim, cell_seed, gate_c_blocked
##
## gate_c_blocked = TRUE means split_residuals=FALSE + tcf="notyettreated":
##   post-Stage-4 code fires Gate C on every fect() call for these cells.
##   The harness records gate_blocked and reports Stage 3 reference coverage.

scenario_grid <- data.frame(
  cell_id          = c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8"),
  cell_index       = 1:8,
  method           = c("ife", "ife", "ife", "ife", "cfe", "cfe", "ife", "ife"),
  tcf              = c("notyettreated", "notyettreated",
                       "notyettreated", "notyettreated",
                       "notyettreated", "notyettreated",
                       "nevertreated",  "nevertreated"),
  dgp_type         = c("ife", "ife", "ife", "ife", "cfe", "cfe", "ife", "ife"),
  N                = c(50L, 50L, 100L, 100L, 50L, 50L, 50L, 50L),
  TT               = rep(20L, 8L),
  Ntr              = c(10L, 10L, 20L, 20L, 10L, 10L, 10L, 10L),
  tau              = rep(1.0, 8L),
  sigma_eps        = rep(1.0, 8L),
  split_residuals  = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE),
  N_sim            = c(1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 500L, 500L),
  cell_seed        = sapply(1:8, derive_cell_seed),
  gate_c_blocked   = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

# Apply N_sim override if provided
if (!is.na(n_sim_override)) {
  n_sim_val <- as.integer(n_sim_override)
  message(sprintf("--n-sim override: all cells set to N_sim=%d", n_sim_val))
  scenario_grid$N_sim <- n_sim_val
}

message(sprintf("Scenario grid: %d cells (D1-D8)", nrow(scenario_grid)))
message(sprintf("Active cells (no Gate C block): %s",
                paste(scenario_grid$cell_id[!scenario_grid$gate_c_blocked],
                      collapse = ", ")))
message(sprintf("Gate C cells (Stage 3 ref used): %s",
                paste(scenario_grid$cell_id[scenario_grid$gate_c_blocked],
                      collapse = ", ")))
total_reps_active <- sum(scenario_grid$N_sim[!scenario_grid$gate_c_blocked])
message(sprintf("Total active fect() calls: %d", total_reps_active))
message(sprintf("Outer mclapply cores: %d", n_cores))

## ---- Single-replication worker for active cells ------------------------------
## This function runs in a forked child process (mclapply). It must:
##   1. Set the seed from the per-rep seed FIRST.
##   2. Generate DGP data with seed=NULL (seed already set above).
##   3. Call fect() with split_residuals per cell spec.
##   4. Return a standardized list from extract_result_debias().
run_one_rep_debias <- function(s, cell_params) {
  rep_seed <- derive_rep_seed(cell_params$cell_seed, s)
  set.seed(rep_seed)

  N         <- cell_params$N
  TT        <- cell_params$TT
  Ntr       <- cell_params$Ntr
  tau       <- cell_params$tau
  sigma_eps <- cell_params$sigma_eps
  dgp_type  <- cell_params$dgp_type
  method    <- cell_params$method
  tcf       <- cell_params$tcf
  split_res <- cell_params$split_residuals

  # --- Generate DGP data -------------------------------------------------------
  simdf <- tryCatch({
    if (dgp_type == "ife") {
      generate_ife_dgp(N = N, TT = TT, Ntr = Ntr, tau = tau,
                       sigma_eps = sigma_eps, r = 2L, seed = NULL)
    } else if (dgp_type == "cfe") {
      generate_cfe_dgp(N = N, TT = TT, Ntr = Ntr, tau = tau,
                       sigma_eps = sigma_eps, r = 2L, beta_z = 1.5,
                       K_levels = 4L, seed = NULL)
    } else {
      stop(paste("Unknown dgp_type:", dgp_type))
    }
  }, error = function(e) {
    structure(paste("DGP error:", conditionMessage(e)),
              class = "try-error", condition = e)
  })

  if (inherits(simdf, "try-error")) {
    return(list(
      covered      = NA_integer_,
      att_est      = NA_real_,
      ci_lower     = NA_real_,
      ci_upper     = NA_real_,
      se_est       = NA_real_,
      failed       = TRUE,
      gate_blocked = FALSE,
      fail_msg     = as.character(simdf)
    ))
  }

  # --- Build fect() index and formula -----------------------------------------
  # IFE cells: index = c("id","time"), formula Y ~ D
  # CFE cells: index = c("id","time","Z"), formula Y ~ D + Z
  use_Z     <- (dgp_type == "cfe")
  index_arg <- if (use_Z) c("id", "time", "Z") else c("id", "time")
  fml       <- if (use_Z) Y ~ D + Z else Y ~ D

  # --- Re-set seed right before fect() for bootstrap reproducibility ----------
  set.seed(rep_seed)

  # --- Call fect() -------------------------------------------------------------
  out <- tryCatch({
    if (tcf == "nevertreated") {
      suppressWarnings(suppressMessages(fect(
        fml,
        data            = simdf,
        index           = index_arg,
        method          = method,
        time.component.from = "nevertreated",
        r               = 2L,
        CV              = FALSE,
        force           = "two-way",
        se              = TRUE,
        vartype         = "parametric",
        nboots          = 200L,
        parallel        = FALSE,
        split_residuals = split_res
      )))
    } else {
      # tcf = "notyettreated" — do NOT pass time.component.from (use default)
      # split_residuals=TRUE bypasses Gate C per post-Stage-4 code
      suppressWarnings(suppressMessages(fect(
        fml,
        data            = simdf,
        index           = index_arg,
        method          = method,
        r               = 2L,
        CV              = FALSE,
        force           = "two-way",
        se              = TRUE,
        vartype         = "parametric",
        nboots          = 200L,
        parallel        = FALSE,
        split_residuals = split_res
      )))
    }
  }, error = function(e) {
    structure(conditionMessage(e), class = "try-error", condition = e)
  })

  extract_result_debias(out, tau_true = tau)
}

## ---- Main simulation loop ----------------------------------------------------
message(sprintf("\nStarting debiasing simulation: %d cells, cores=%d",
                nrow(scenario_grid), n_cores))
message(sprintf("Master seed: %d\n", MASTER_SEED))

wall_start <- proc.time()["elapsed"]

all_raw_results   <- vector("list", nrow(scenario_grid))
all_summary_rows  <- vector("list", nrow(scenario_grid))

for (gi in seq_len(nrow(scenario_grid))) {
  cell    <- as.list(scenario_grid[gi, ])
  cell_id <- cell$cell_id
  N_sim   <- cell$N_sim

  # --- --only-cells filter ---------------------------------------------------
  if (!is.null(only_cells) && !(cell_id %in% only_cells)) {
    message(sprintf("[%d/8] SKIP (not in --only-cells): %s", gi, cell_id))
    next
  }

  # --- Load existing checkpoint -----------------------------------------------
  chk_path <- cell_rds_path(cell_id, out_dir)
  if (file.exists(chk_path)) {
    message(sprintf("[%d/8] SKIP (checkpoint exists): %s -> %s",
                    gi, cell_id, chk_path))
    saved <- readRDS(chk_path)
    all_raw_results[[gi]]  <- saved$raw
    all_summary_rows[[gi]] <- saved$summary_row
    next
  }

  # --- Gate C blocked cells (D1, D3, D5) -------------------------------------
  if (isTRUE(cell$gate_c_blocked)) {
    message(sprintf("[%d/8] GATE-C BLOCKED: %s (split_residuals=FALSE, tcf=notyettreated)",
                    gi, cell_id))

    stage3_ref <- STAGE3_REF[[cell_id]]

    # Run a small probe (R=3) to confirm Gate C fires; then set gate_blocked
    probe_reps <- parallel::mclapply(seq_len(3L), function(s) {
      set.seed(derive_rep_seed(cell$cell_seed, s))
      simdf <- tryCatch({
        generate_ife_dgp(N = cell$N, TT = cell$TT, Ntr = cell$Ntr,
                         tau = cell$tau, sigma_eps = cell$sigma_eps,
                         r = 2L, seed = NULL)
      }, error = function(e)
        structure(paste("DGP error:", conditionMessage(e)),
                  class = "try-error"))
      if (inherits(simdf, "try-error")) {
        return(list(failed = TRUE, gate_blocked = FALSE,
                    fail_msg = as.character(simdf)))
      }
      use_Z     <- (cell$dgp_type == "cfe")
      index_arg <- if (use_Z) c("id", "time", "Z") else c("id", "time")
      fml       <- if (use_Z) Y ~ D + Z else Y ~ D
      out_probe <- tryCatch(
        suppressWarnings(suppressMessages(fect(
          fml,
          data            = simdf,
          index           = index_arg,
          method          = cell$method,
          r               = 2L,
          CV              = FALSE,
          force           = "two-way",
          se              = TRUE,
          vartype         = "parametric",
          nboots          = 5L,
          parallel        = FALSE,
          split_residuals = FALSE
        ))),
        error = function(e)
          structure(conditionMessage(e), class = "try-error")
      )
      extract_result_debias(out_probe, tau_true = cell$tau)
    }, mc.cores = min(3L, n_cores), mc.set.seed = FALSE)

    n_gate_probe <- sum(sapply(probe_reps,
                               function(r) isTRUE(r$gate_blocked)))
    confirmed_blocked <- (n_gate_probe == 3L)
    message(sprintf("  Gate C probe: %d/3 reps blocked. confirmed=%s",
                    n_gate_probe, confirmed_blocked))

    # Build a synthetic summary row using Stage 3 reference
    agg_blocked <- list(
      n_total   = 0L,
      n_valid   = 0L,
      n_failed  = 0L,
      n_gate    = 0L,
      cov_rate  = stage3_ref$coverage_rate,
      mc_se_cov = if (!is.na(stage3_ref$coverage_rate) &&
                       stage3_ref$N_sim > 0L) {
        sqrt(stage3_ref$coverage_rate *
             (1.0 - stage3_ref$coverage_rate) / stage3_ref$N_sim)
      } else {
        NA_real_
      },
      mean_att  = NA_real_,
      bias      = NA_real_,
      sd_att    = NA_real_,
      mean_se   = NA_real_,
      se_ratio  = NA_real_,
      fail_rate = NA_real_,
      gate_rate = NA_real_
    )

    summary_row_blocked <- data.frame(
      cell_id         = cell_id,
      cell_index      = cell$cell_index,
      method          = cell$method,
      tcf             = cell$tcf,
      dgp_type        = cell$dgp_type,
      N               = cell$N,
      TT              = cell$TT,
      Ntr             = cell$Ntr,
      tau             = cell$tau,
      sigma_eps       = cell$sigma_eps,
      split_residuals = cell$split_residuals,
      N_sim           = N_sim,
      N_sim_valid     = 0L,
      N_sim_failed    = 0L,
      N_gate_blocked  = 0L,
      coverage        = round(agg_blocked$cov_rate, 4L),
      mc_se_cov       = round(agg_blocked$mc_se_cov, 5L),
      mean_att        = NA_real_,
      bias            = NA_real_,
      sd_att          = NA_real_,
      mean_se         = NA_real_,
      se_ratio        = NA_real_,
      fail_pct        = NA_real_,
      gate_pct        = NA_real_,
      elapsed_sec     = 0.0,
      gate_c_blocked  = TRUE,
      stage3_ref      = stage3_ref$coverage_rate,
      stage3_N_sim    = stage3_ref$N_sim,
      stage3_note     = stage3_ref$note,
      stringsAsFactors = FALSE
    )

    save_cell_checkpoint(
      list(raw = list(cell = cell, reps = probe_reps,
                      agg = agg_blocked, elapsed = 0.0,
                      gate_confirmed = confirmed_blocked,
                      stage3_ref = stage3_ref),
           summary_row = summary_row_blocked),
      cell_id = cell_id,
      out_dir = out_dir
    )
    message(sprintf("  Checkpoint saved: cell-%s.rds (Gate C, Stage 3 ref=%.3f)",
                    cell_id,
                    ifelse(is.na(stage3_ref$coverage_rate), NA,
                           stage3_ref$coverage_rate)))

    all_raw_results[[gi]]  <- list(cell = cell, reps = probe_reps,
                                   agg = agg_blocked, elapsed = 0.0,
                                   gate_confirmed = confirmed_blocked,
                                   stage3_ref = stage3_ref)
    all_summary_rows[[gi]] <- summary_row_blocked
    next
  }

  # --- Active cells (D2, D4, D6, D7, D8) ------------------------------------
  cell_label <- sprintf(
    "[%d/8] %s | method=%s tcf=%s dgp=%s N=%d Ntr=%d split=%s N_sim=%d",
    gi, cell_id, cell$method, cell$tcf, cell$dgp_type,
    cell$N, cell$Ntr,
    as.character(cell$split_residuals),
    N_sim
  )
  message(cell_label)

  t_cell_start <- proc.time()["elapsed"]

  reps <- parallel::mclapply(
    seq_len(N_sim),
    FUN          = run_one_rep_debias,
    cell_params  = cell,
    mc.cores     = n_cores,
    mc.set.seed  = FALSE   # seeds managed explicitly inside run_one_rep_debias
  )

  t_cell_end <- proc.time()["elapsed"]
  cell_elapsed <- t_cell_end - t_cell_start

  # Aggregate
  agg <- aggregate_cell_debias(reps, tau_true = cell$tau)

  # Progress message
  fail_flag <- if (!is.na(agg$fail_rate) && agg$fail_rate > 0.05)
    "  *** NOTABLE FAILURE RATE ***" else ""
  message(sprintf(
    "  Done: coverage=%.4f (MC_SE=%.4f), bias=%.5f, SE_ratio=%.4f, fail=%.1f%%%s | %.1f sec",
    agg$cov_rate, agg$mc_se_cov,
    agg$bias, agg$se_ratio,
    100.0 * agg$fail_rate,
    fail_flag,
    cell_elapsed
  ))

  summary_row <- data.frame(
    cell_id         = cell_id,
    cell_index      = cell$cell_index,
    method          = cell$method,
    tcf             = cell$tcf,
    dgp_type        = cell$dgp_type,
    N               = cell$N,
    TT              = cell$TT,
    Ntr             = cell$Ntr,
    tau             = cell$tau,
    sigma_eps       = cell$sigma_eps,
    split_residuals = cell$split_residuals,
    N_sim           = N_sim,
    N_sim_valid     = agg$n_valid,
    N_sim_failed    = agg$n_failed,
    N_gate_blocked  = agg$n_gate,
    coverage        = round(agg$cov_rate, 4L),
    mc_se_cov       = round(agg$mc_se_cov, 5L),
    mean_att        = round(agg$mean_att, 5L),
    bias            = round(agg$bias, 5L),
    sd_att          = round(agg$sd_att, 5L),
    mean_se         = round(agg$mean_se, 5L),
    se_ratio        = round(agg$se_ratio, 4L),
    fail_pct        = round(100.0 * agg$fail_rate, 2L),
    gate_pct        = round(100.0 * agg$gate_rate, 2L),
    elapsed_sec     = round(cell_elapsed, 1L),
    gate_c_blocked  = FALSE,
    stage3_ref      = NA_real_,
    stage3_N_sim    = NA_integer_,
    stage3_note     = NA_character_,
    stringsAsFactors = FALSE
  )

  raw_result <- list(
    cell    = cell,
    reps    = reps,
    agg     = agg,
    elapsed = cell_elapsed
  )

  save_cell_checkpoint(
    list(raw = raw_result, summary_row = summary_row),
    cell_id = cell_id,
    out_dir = out_dir
  )
  message(sprintf("  Checkpoint saved: cell-%s.rds", cell_id))

  all_raw_results[[gi]]  <- raw_result
  all_summary_rows[[gi]] <- summary_row

  # ETA projection
  n_done <- sum(!sapply(all_raw_results, is.null))
  if (n_done >= 2L) {
    elapsed_so_far  <- proc.time()["elapsed"] - wall_start
    avg_per_cell    <- elapsed_so_far / n_done
    n_remaining     <- sum(sapply(seq_len(nrow(scenario_grid)), function(j) {
      cid_j    <- scenario_grid$cell_id[j]
      in_filt  <- is.null(only_cells) || (cid_j %in% only_cells)
      not_done <- !file.exists(cell_rds_path(cid_j, out_dir))
      in_filt && not_done && !scenario_grid$gate_c_blocked[j]
    }))
    eta_sec <- avg_per_cell * n_remaining
    message(sprintf("  ETA: ~%.0f min remaining (%.1f min elapsed)",
                    eta_sec / 60.0, elapsed_so_far / 60.0))
  }

  flush.console()
}

## ---- Load any checkpointed cells not in memory (partial restart support) -----
message("\nAggregating all available cell checkpoints...")
for (gi in seq_len(nrow(scenario_grid))) {
  if (!is.null(all_raw_results[[gi]])) next
  cid_j    <- scenario_grid$cell_id[gi]
  chk_path <- cell_rds_path(cid_j, out_dir)
  if (file.exists(chk_path)) {
    saved <- readRDS(chk_path)
    all_raw_results[[gi]]  <- saved$raw
    all_summary_rows[[gi]] <- saved$summary_row
    message(sprintf("  Loaded checkpoint: %s", cid_j))
  }
}

wall_end      <- proc.time()["elapsed"]
total_elapsed <- wall_end - wall_start

## ---- Assemble final summary data.frame ---------------------------------------
valid_summary <- Filter(Negate(is.null), all_summary_rows)
valid_raw     <- Filter(Negate(is.null), all_raw_results)

if (length(valid_summary) == 0L) {
  message("No completed cells — nothing to aggregate. Exiting.")
  quit(status = 0L)
}

summary_df <- do.call(rbind, valid_summary)
rownames(summary_df) <- NULL

## ---- Save results -----------------------------------------------------------
results_obj <- list(
  scenario_grid = scenario_grid,
  summary       = summary_df,
  raw_results   = valid_raw,
  stage3_ref    = STAGE3_REF,
  meta = list(
    master_seed       = MASTER_SEED,
    n_cores           = n_cores,
    fect_path         = fect_path,
    fect_git_commit   = paste(fect_git_commit, collapse = " "),
    r_version         = R.version.string,
    total_elapsed_sec = total_elapsed,
    run_date          = Sys.time()
  )
)

rds_path <- file.path(out_dir, "debias-results.rds")
saveRDS(results_obj, rds_path)
message(sprintf("\nResults saved to: %s", rds_path))

summary_rds <- file.path(out_dir, "debias-summary.rds")
saveRDS(summary_df, summary_rds)
message(sprintf("Summary saved to: %s", summary_rds))

## ---- Go/No-Go determination --------------------------------------------------
## Extract key cell coverage values for the decision rule
get_cov <- function(cid) {
  row <- summary_df[summary_df$cell_id == cid, ]
  if (nrow(row) == 0L) return(NA_real_)
  row$coverage[1L]
}

D2_cov <- get_cov("D2")
D4_cov <- get_cov("D4")
D8_cov <- get_cov("D8")

go_cond_1  <- !is.na(D2_cov) && D2_cov >= 0.92
go_cond_2  <- !is.na(D4_cov) && D4_cov >= 0.92
go_cond_3  <- !is.na(D8_cov) && D8_cov >= 0.92 && D8_cov <= 0.98
partial    <- !go_cond_1 && go_cond_2 && go_cond_3

verdict <- if (go_cond_1 && go_cond_2 && go_cond_3) {
  "GO"
} else if (partial) {
  "PARTIAL"
} else {
  "NO-GO"
}

## ---- Final console output ----------------------------------------------------
message("\n\n========== DEBIASING COVERAGE SIMULATION RESULTS ==========\n")
message(sprintf("Total wall time:  %.1f minutes", total_elapsed / 60.0))
message(sprintf("R:                %s", R.version.string))
message(sprintf("fect HEAD:        %s", paste(fect_git_commit, collapse = " ")))
message(sprintf("Cores used:       %d", n_cores))
message(sprintf("Master seed:      %d\n", MASTER_SEED))

message("=== COVERAGE TABLE (all 8 cells) ===")
header <- sprintf("%-4s  %-6s  %-14s  %-4s  %5s  %5s  %-7s  %9s  %7s  %8s  %9s  %8s  %7s  %5s  %5s",
                  "Cell", "method", "tcf", "dgp", "N", "Nsim",
                  "split", "coverage", "MC_SE",
                  "bias", "emp_SD", "mean_SE", "SE_ratio",
                  "fail%", "gate%")
message(header)
message(paste(rep("-", nchar(header)), collapse = ""))

for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cov_str <- if (is.na(r$coverage)) {
    sprintf("  (Stage3: %s)",
            if (!is.na(r$stage3_ref))
              sprintf("%.3f [n=%d]", r$stage3_ref, r$stage3_N_sim)
            else "NA")
  } else {
    sprintf("%9.4f", r$coverage)
  }
  message(sprintf("%-4s  %-6s  %-14s  %-4s  %5d  %5d  %-7s  %s  %7.5f  %8.5f  %9.5f  %8.5f  %7.4f  %5.1f  %5.1f",
                  r$cell_id, r$method, r$tcf, r$dgp_type, r$N, r$N_sim,
                  as.character(r$split_residuals),
                  cov_str,
                  ifelse(is.na(r$mc_se_cov), 0, r$mc_se_cov),
                  ifelse(is.na(r$bias),      0, r$bias),
                  ifelse(is.na(r$sd_att),    0, r$sd_att),
                  ifelse(is.na(r$mean_se),   0, r$mean_se),
                  ifelse(is.na(r$se_ratio),  0, r$se_ratio),
                  ifelse(is.na(r$fail_pct),  0, r$fail_pct),
                  ifelse(is.na(r$gate_pct),  0, r$gate_pct)))
}

message("\n=== PAIRED COMPARISON (baseline vs debiased) ===")
pairs <- list(c("D1", "D2"), c("D3", "D4"), c("D5", "D6"), c("D7", "D8"))
for (p in pairs) {
  r_base  <- summary_df[summary_df$cell_id == p[1], ]
  r_split <- summary_df[summary_df$cell_id == p[2], ]
  base_cov  <- if (nrow(r_base) > 0L) {
    if (!is.na(r_base$coverage[1L])) r_base$coverage[1L] else r_base$stage3_ref[1L]
  } else NA_real_
  split_cov <- if (nrow(r_split) > 0L) r_split$coverage[1L] else NA_real_
  delta     <- if (!is.na(base_cov) && !is.na(split_cov)) split_cov - base_cov else NA_real_
  base_label  <- if (nrow(r_base) > 0L && isTRUE(r_base$gate_c_blocked[1L]))
    "(Stage3 ref)" else "(new)"
  message(sprintf("  %s vs %s:  baseline=%.4f %s  debiased=%.4f  delta=%+.4f",
                  p[1], p[2],
                  ifelse(is.na(base_cov),  NA_real_, base_cov),  base_label,
                  ifelse(is.na(split_cov), NA_real_, split_cov),
                  ifelse(is.na(delta),     NA_real_, delta)))
}

message("\n=== GO/NO-GO VERDICT ===")
message(sprintf("  D2 (ife+nt N=50  split=TRUE):  coverage=%.4f  (threshold >= 0.92)  %s",
                ifelse(is.na(D2_cov), NA_real_, D2_cov),
                if (go_cond_1) "PASS" else "FAIL"))
message(sprintf("  D4 (ife+nt N=100 split=TRUE):  coverage=%.4f  (threshold >= 0.92)  %s",
                ifelse(is.na(D4_cov), NA_real_, D4_cov),
                if (go_cond_2) "PASS" else "FAIL"))
message(sprintf("  D8 (ife+nev N=50 split=TRUE):  coverage=%.4f  (threshold [0.92,0.98])  %s",
                ifelse(is.na(D8_cov), NA_real_, D8_cov),
                if (go_cond_3) "PASS" else "FAIL"))
message(sprintf("\n  VERDICT: %s", verdict))

message("\n=== SEED DOCUMENTATION ===")
message(sprintf("  Master seed:      %d", MASTER_SEED))
message("  Cell seed formula: MASTER_SEED + (cell_index - 1) * 10000")
for (i in seq_len(nrow(scenario_grid))) {
  message(sprintf("    %s (index=%d): seed=%d",
                  scenario_grid$cell_id[i],
                  scenario_grid$cell_index[i],
                  scenario_grid$cell_seed[i]))
}
message("  Rep seed formula:  cell_seed + rep_index")
message("  RNG type:          Mersenne-Twister (R default, no RNGkind() call)")

message(sprintf("\nSimulation complete. All results in: %s", out_dir))
