## Addendum to ar1-reference-gen.R (round 3, NOT run by the tests): the prototype's treated_rank() on
## gsynth fits of P4 with max.iteration = 3 (test N5, "plain path untouched"). Adds
## rules$P4_maxit3 to ar1-reference.rds. Same prototype checkout and requirements as ar1-reference-gen.R
## (xuyiqing/gsynth-revisit wip/template-v2-2026-10-16, e1eebae; R/rank_rule.R), fits from the base build
## (6c9c16c, library on R_LIBS), calls: method = "gsynth", force = "two-way", CV = FALSE, r = 0..5,
## tol = 1e-10, max.iteration = 3, parallel = FALSE (fect's own "EM did not converge" warnings suppressed).
##   PROTO_DIR=<checkout> R_LIBS=<base lib> Rscript ar1-reference-gen-round3.R
## Generated 2026-10-16 (tester round 3).
suppressPackageStartupMessages(library(fect))
gen_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
proto <- Sys.getenv("PROTO_DIR"); stopifnot(nzchar(proto))
source(file.path(gen_dir, "..", "helper-ar1.R"))
source(file.path(proto, "R", "rank_rule.R"))
fl <- lapply(0:5, function(r) suppressWarnings(ar1_fect(ar1_P4(), CV = FALSE, r = r, max.iteration = 3)))
names(fl) <- 0:5
ref <- readRDS(file.path(gen_dir, "ar1-reference.rds"))
ref$rules$P4_maxit3 <- treated_rank(fl)
ref$meta$round3 <- list(note = "rules$P4_maxit3 added", fect_for_fits = as.character(packageVersion("fect")), date = "2026-10-16")
saveRDS(ref, file.path(gen_dir, "ar1-reference.rds"), version = 2)
print(ref$rules$P4_maxit3$table); print(ref$rules$P4_maxit3$r_loo)
