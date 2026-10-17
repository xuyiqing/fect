## C3: cv.method = "pretreatment" reproduces the template's treated-unit-rule picks on all 49 cases.
##
## Reference (test-spec C3, amended 2026-10-16): <GSYNTH_REVISIT_DATA>/synthesis/runs/2026-10-16-rank-figure/rank_estimates.csv
## (columns paper, case, r_rule, r_fwd, rule_note; written by the prototype's scripts/rank_estimates.R at ae489b6) and the per-rank
## scores papers/<paper>/output/<case>/tables/rank-rule.csv (columns r, loo, fwd), computed by the prototype from cached gsynth point fits.
##
## Heavy and needs data that is not in the package: skip_on_cran(), and it skips unless these environment variables are set
##   GSYNTH_REVISIT_DATA  folder (or a local copy of it) with synthesis/ and papers/<paper>/ (spec, data, output/<case>/tables)
##   GSYNTH_REVISIT_PROTO a checkout/archive of the prototype at ae489b6 (R/paths.R, R/template.R, papers/<paper>/00_spec.R)
## Optional: AR1_C3_PAPERS (comma-separated paper ids to restrict to), AR1_C3_LOG (csv to append one row per case).
## Eibl cases take minutes; run with  Rscript -e 'testthat::test_file("tests/testthat/test-rank-cases.R")'.

skip_on_cran()
data_root <- Sys.getenv("GSYNTH_REVISIT_DATA"); proto <- Sys.getenv("GSYNTH_REVISIT_PROTO")
skip_if(!nzchar(data_root) || !nzchar(proto) || !dir.exists(data_root) || !dir.exists(proto),
        "set GSYNTH_REVISIT_DATA and GSYNTH_REVISIT_PROTO (prototype ae489b6) to run the 49-case rank check")
est_file <- file.path(data_root, "synthesis", "runs", "2026-10-16-rank-figure", "rank_estimates.csv")
skip_if(!file.exists(est_file), "rank_estimates.csv not found")

est <- read.csv(est_file, stringsAsFactors = FALSE)
est <- est[est$status == "ok", ]
only <- strsplit(Sys.getenv("AR1_C3_PAPERS"), ",")[[1]]
if (length(only)) est <- est[est$paper %in% only, ]

note_ints <- function(note) {      # counts and percentages the note reports
  m <- regmatches(note, regexec("ranks 0-([0-9]+); leave one out scores ([0-9]+) of ([0-9]+) treated units \\(([0-9]+)% of the ATT weight\\), forward ([0-9]+) \\(([0-9]+)%\\)", note))[[1]]
  setNames(as.integer(m[-1]), c("k", "n_loo", "n_tr", "p_loo", "n_fwd", "p_fwd"))
}

load_env <- local({
  e <- NULL
  function() {
    if (is.null(e)) {
      Sys.setenv(GSYNTH_REVISIT_DATA = data_root)
      e <<- globalenv()          # the case specs call paper_path() and friends from the global environment
      withr::with_dir(proto, {
        sys.source("R/paths.R", envir = e)
        suppressMessages(sys.source("R/template.R", envir = e))
      })
    }
    e
  }
})

case_args <- function(paper, case, k) {
  e <- load_env()
  withr::with_dir(proto, {
    pe <- new.env(); local(source(file.path("papers", paper, "00_spec.R"), local = TRUE), envir = pe)
    cell <- pe$cells[[which(vapply(pe$cells, `[[`, "", "id") == case)]]
    s <- e$cell_spec(pe$spec, cell)
    args <- list(formula = e$cell_formula(s, cell), data = e$cell_data(s, cell), index = s$index, method = "gsynth",
                 force = s$force, CV = TRUE, r = c(0, k), cv.method = "pretreatment", se = FALSE, parallel = FALSE, seed = s$seed)
    if (!is.null(s$min.T0)) args$min.T0 <- s$min.T0
    c(args, s$fect_args)
  })
}

rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-300)

test_that("C3: reference file lists 49 cases", {
  if (!length(only)) expect_equal(nrow(est), 49)
})

for (i in seq_len(nrow(est))) local({
  row <- est[i, ]; id <- sprintf("%s/%s", row$paper, row$case)
  test_that(sprintf("C3: %s rank rule", id), {
    ni <- note_ints(row$rule_note); k <- ni[["k"]]
    rr <- read.csv(file.path(data_root, "papers", row$paper, "output", row$case, "tables", "rank-rule.csv"))
    expect_identical(rr$r[which.min(rr$loo)], as.integer(row$r_rule), info = paste(id, "reference table argmin = r_rule"))
    t0 <- Sys.time()
    args <- suppressWarnings(suppressMessages(case_args(row$paper, row$case, k)))   # readxl cell-type warnings from the paper's data loader
    ## fect's own "EM did not converge within max.iteration" warning on a case's data is recorded in the log, not re-raised
    em_warn <- FALSE
    fit <- withCallingHandlers(suppressMessages(do.call(fect::fect, args)),
                               warning = function(w) if (grepl("EM did not converge", conditionMessage(w), fixed = TRUE)) {
                                 em_warn <<- TRUE; invokeRestart("muffleWarning") })
    ct <- fit$cv.treated; tab <- ct$table
    ## 1. picks (ties at EM tolerance are reported, not failed)
    loo_ref <- setNames(rr$loo, rr$r)
    tie <- !identical(as.integer(fit$r.cv), as.integer(row$r_rule)) &&
      rel(loo_ref[[as.character(fit$r.cv)]], loo_ref[[as.character(row$r_rule)]]) < 1e-6
    if (tie) message(sprintf("C3 tie at EM tolerance: %s build %d vs reference %d", id, fit$r.cv, row$r_rule))
    if (!tie) {
      expect_identical(as.integer(fit$r.cv), as.integer(row$r_rule), info = id)
      expect_identical(as.integer(ct$r_loo), as.integer(row$r_rule), info = id)
    }
    ## 2. forward pick
    expect_identical(as.integer(ct$r_fwd), as.integer(row$r_fwd), info = paste(id, "r_fwd"))
    ## 3. unit counts and percentages
    expect_identical(as.integer(c(ct$n_tr, ct$n_loo, round(100 * ct$share_loo), ct$n_fwd, round(100 * ct$share_fwd))),
                     as.integer(c(ni[["n_tr"]], ni[["n_loo"]], ni[["p_loo"]], ni[["n_fwd"]], ni[["p_fwd"]])), info = paste(id, "counts"))
    ## 4. scores at every rank 0..k, relative 1e-6
    rr <- rr[rr$r <= k, ]
    expect_identical(sort(tab$r), sort(rr$r[is.finite(rr$loo)]), info = paste(id, "scored ranks"))
    m <- match(rr$r, tab$r)
    dl <- rel(tab$loo[m], rr$loo); dl <- dl[is.finite(dl)]
    fin <- is.finite(rr$fwd) & !is.na(m); df <- rel(tab$fwd[m][fin], rr$fwd[fin]); df <- df[is.finite(df)]
    max_l <- max(c(dl, 0)); max_f <- max(c(df, 0))
    expect_lte(max_l, 1e-6); expect_lte(max_f, 1e-6)
    lg <- Sys.getenv("AR1_C3_LOG")
    if (nzchar(lg)) cat(sprintf("%s,%s,%d,%d,%d,%d,%d,%d,%d,%.3e,%.3e,%s,%s,%.0f\n", row$paper, row$case, k, row$r_rule, fit$r.cv, row$r_fwd,
                                 ct$r_fwd, ct$n_loo, ct$n_tr, max_l, max_f, tie, em_warn, as.numeric(Sys.time() - t0, units = "secs")),
                        file = lg, append = TRUE)
  })
})
