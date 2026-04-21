## make-paraboot-baseline.R
## Pre-refactor baseline fixture capture helper (test-spec.md §1)
##
## Run this script from the fect project root BEFORE applying the refactor:
##   Rscript tests/testthat/make-paraboot-baseline.R
##
## Saves: tests/testthat/fixtures/paraboot-baseline.rds
##
## Seed convention: set.seed(42) for dataset; set.seed(2026) before each fect() call.
## parallel=FALSE is MANDATORY for reproducibility.

library(fect)

make_fixture_data <- function(seed = 42) {
  set.seed(seed)
  N   <- 50
  TT  <- 20
  T0  <- 12
  Ntr <- 10
  id   <- rep(1:N, each = TT)
  time <- rep(1:TT, N)
  D    <- as.integer(id <= Ntr & time >= T0)
  alpha_i <- rep(rnorm(N, 0, 1), each = TT)
  xi_t    <- rep(rnorm(TT, 0, 0.5), N)
  Y <- alpha_i + xi_t + 2 * D + rnorm(N * TT, 0, 0.5)
  data.frame(id = id, time = time, Y = Y, D = D)
}

d <- make_fixture_data(seed = 42)

## Fixture calls — test-spec.md §1.2
set.seed(2026)
gsynth_para_result <- suppressWarnings(suppressMessages(fect(
  Y ~ D, data=d, index=c("id","time"),
  method="gsynth", r=1, se=TRUE, vartype="parametric", nboots=50, CV=FALSE, parallel=FALSE
)))

set.seed(2026)
ife_nt_para_result <- suppressWarnings(suppressMessages(fect(
  Y ~ D, data=d, index=c("id","time"),
  method="ife", r=1, se=TRUE, vartype="parametric", nboots=50, CV=FALSE, parallel=FALSE
)))

set.seed(2026)
ife_nev_para_result <- suppressWarnings(suppressMessages(fect(
  Y ~ D, data=d, index=c("id","time"),
  method="ife", time.component.from="nevertreated", r=1, se=TRUE,
  vartype="parametric", nboots=50, CV=FALSE, parallel=FALSE
)))

set.seed(2026)
cfe_nev_para_result <- suppressWarnings(suppressMessages(fect(
  Y ~ D, data=d, index=c("id","time"),
  method="cfe", time.component.from="nevertreated", r=1, se=TRUE,
  vartype="parametric", nboots=50, CV=FALSE, parallel=FALSE
)))

set.seed(2026)
cfe_nt_para_result <- suppressWarnings(suppressMessages(fect(
  Y ~ D, data=d, index=c("id","time"),
  method="cfe", r=0, se=TRUE, vartype="parametric", nboots=50, CV=FALSE, parallel=FALSE
)))

set.seed(2026)
ife_nt_boot_result <- suppressWarnings(suppressMessages(fect(
  Y ~ D, data=d, index=c("id","time"),
  method="ife", r=1, se=TRUE, vartype="bootstrap", nboots=50, CV=FALSE, parallel=FALSE
)))

set.seed(2026)
ife_nev_boot_result <- suppressWarnings(suppressMessages(fect(
  Y ~ D, data=d, index=c("id","time"),
  method="ife", time.component.from="nevertreated", r=1, se=TRUE,
  vartype="bootstrap", nboots=50, CV=FALSE, parallel=FALSE
)))

set.seed(2026)
gsynth_boot_result <- suppressWarnings(suppressMessages(fect(
  Y ~ D, data=d, index=c("id","time"),
  method="gsynth", r=1, se=TRUE, vartype="bootstrap", nboots=50, CV=FALSE, parallel=FALSE
)))

set.seed(2026)
ife_jk_result <- suppressWarnings(suppressMessages(fect(
  Y ~ D, data=d, index=c("id","time"),
  method="ife", r=1, se=TRUE, vartype="jackknife", CV=FALSE, parallel=FALSE
)))

set.seed(2026)
gsynth_jk_result <- suppressWarnings(suppressMessages(fect(
  Y ~ D, data=d, index=c("id","time"),
  method="gsynth", r=1, se=TRUE, vartype="jackknife", CV=FALSE, parallel=FALSE
)))

## Save
baseline <- list(
  gsynth_para  = gsynth_para_result,
  ife_nt_para  = ife_nt_para_result,
  ife_nev_para = ife_nev_para_result,
  cfe_nev_para = cfe_nev_para_result,
  cfe_nt_para  = cfe_nt_para_result,
  ife_nt_boot  = ife_nt_boot_result,
  ife_nev_boot = ife_nev_boot_result,
  gsynth_boot  = gsynth_boot_result,
  ife_jk       = ife_jk_result,
  gsynth_jk    = gsynth_jk_result
)

dir.create("tests/testthat/fixtures", showWarnings = FALSE, recursive = TRUE)
saveRDS(baseline, "tests/testthat/fixtures/paraboot-baseline.rds")
message("Fixture saved to tests/testthat/fixtures/paraboot-baseline.rds")
