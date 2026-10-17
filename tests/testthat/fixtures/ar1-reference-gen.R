## Reference generator for test-cv-treated.R, test-ar1-joint.R (NOT run by the tests).
##
## Needs the gsynth-revisit prototype checkout (xuyiqing/gsynth-revisit, branch
## wip/template-v2-2026-10-16, commit e1eebae) and a fect that is NOT the build
## under test (the base 6c9c16c is used: the gsynth fits with CV = FALSE that
## feed the prototype's treated_rank() are unchanged by the AR(1) work).
##
##   PROTO_DIR=~/GitHub/gsynth-revisit-checkout  Rscript ar1-reference-gen.R
##
## where PROTO_DIR holds R/rank_rule.R, sim/ar1_options.R, sim/motivating_pretrend.R,
## sim/rank_rules.R (extract with `git show origin/wip/template-v2-2026-10-16:<path>`).
## Output: ar1-reference.rds beside this file. Every pinned number in the tests
## comes from here (prototype functions run unchanged). Generated 2026-10-16.
##
## Contents
##   P1$Y            20 x 35 outcome matrix (make_design(30, 5, 15, 5, 0.5, seed = 7); set.seed(1);
##                   ar1_paths(35, 20, 0.6)); sum(Y) = 361.4431430351
##   P1$levels_att   levels_att(Y, co, tr, 15, r), r = 0..3
##   P1$filtered     filtered_att(Y, co, tr, 15, r, x), r = 0..3 (rows) by x = 0, 0.3, 0.6, 1
##   P1$x_joint      joint x at r = 0..3 (grid 0.05 on [0, 1], optimize on the bracket, better of the two)
##   P1$path         the alternation: round 0 = levels rule's rank; x_k at r_{k-1}; r_k by the rule
##                   on the filtered panel with T0 - 1 pretreatment periods; stop at a repeated rank
##   P1$rank_scores  the simulation's rank_scores() (loo_cell, fwd_cell rows) on levels, ranks 0..5
##   rules$P1..P5    treated_rank() of the prototype on gsynth fits at r = 0..5 (CV = FALSE)
##   boundary        joint x on white-noise and random-walk panels (prototype ssr_filtered)
##   r10             the prototype's leave-one-out pick (rank_scores) on ar1_mot_panel(seed), seeds 1..40
suppressPackageStartupMessages(library(fect))
gen_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
proto <- Sys.getenv("PROTO_DIR"); stopifnot(nzchar(proto))
source(file.path(gen_dir, "..", "helper-ar1.R"))
source(file.path(proto, "sim", "ar1_options.R"))
source(file.path(proto, "sim", "motivating_pretrend.R")); source(file.path(proto, "sim", "rank_rules.R"))
source(file.path(proto, "R", "rank_rule.R"))
## ar1_options.R and motivating_example.R both define make_design(); re-source the AR(1) one last
source(file.path(proto, "sim", "ar1_options.R"))

## --- P1 -----------------------------------------------------------------------------------
des <- make_design(30, 5, 15, 5, 0.5, seed = 7)
set.seed(1); E <- ar1_paths(35, des$T, 0.6); Y <- des$mu0 + E
co <- 1:30; tr <- 31:35; T0 <- 15
stopifnot(abs(sum(Y) - 361.4431430351) < 1e-9)
P1 <- list(Y = Y, co = co, tr = tr, T0 = T0)
P1$levels_att <- vapply(0:3, function(r) levels_att(Y, co, tr, T0, r), 0)
P1$filtered <- t(vapply(0:3, function(r) vapply(c(0, 0.3, 0.6, 1), function(x) filtered_att(Y, co, tr, T0, r, x), 0), numeric(4)))
dimnames(P1$filtered) <- list(r = 0:3, x = c(0, 0.3, 0.6, 1))

## joint x: the E3 procedure
x_ref <- function(Yc, r) {
  g <- seq(0, 1, by = 0.05); s <- vapply(g, function(x) ssr_filtered(Yc, x, r), 0); k <- which.min(s)
  f <- function(x) ssr_filtered(Yc, x, r)
  xo <- optimize(f, c(g[max(1, k - 1)], g[min(21, k + 1)]), tol = 1e-5)$minimum
  if (f(xo) <= f(g[k])) xo else g[k]
}
P1$x_joint <- vapply(0:3, function(r) x_ref(Y[, co], r), 0)
P1$att_at_joint <- vapply(0:3, function(r) filtered_att(Y, co, tr, T0, r, P1$x_joint[r + 1]), 0)   # filtered_att at the joint x, r = 0..3
rule_pick <- function(Yp, T0p) unname(pick_row(rank_scores(Yp, co, tr, T0p, 5))["loo_cell"])
P1$rank_scores <- rank_scores(Y, co, tr, T0, 5)
## alternation (cap 5 rounds)
r_prev <- rule_pick(Y, T0); path <- data.frame(round = 0L, x = NA_real_, r = r_prev, att = NA_real_)
for (k in 1:5) {
  xk <- x_ref(Y[, co], r_prev)
  V <- Y[-1, , drop = FALSE] - xk * Y[-nrow(Y), , drop = FALSE]
  rk <- rule_pick(V, T0 - 1)
  path <- rbind(path, data.frame(round = k, x = xk, r = rk, att = filtered_att(Y, co, tr, T0, rk, xk)))
  if (rk == r_prev) break
  r_prev <- rk
}
P1$path <- path
## one-round pipeline: x at the levels rule's rank, rank on the filtered panel, estimate
P1$one_round <- path[path$round == 1, ]

## --- the rule on the other test panels (gsynth fits, CV = FALSE) ------------------------------
fits_for <- function(d, formula = y ~ d, min.T0 = NULL) {
  fl <- lapply(0:5, function(r) ar1_fect(d, formula, CV = FALSE, r = r, min.T0 = min.T0))
  names(fl) <- 0:5; fl
}
long1 <- data.frame(id = rep(1:35, each = 20), time = rep(1:20, 35), y = as.vector(Y),
                    d = as.integer(rep(1:35, each = 20) > 30 & rep(1:20, 35) > 15))
rules <- list(P1 = treated_rank(fits_for(long1)),
              P2 = treated_rank(fits_for(ar1_P2())),
              P3 = treated_rank(fits_for(ar1_P3())),
              P4 = treated_rank(fits_for(ar1_P4())),
              P5 = treated_rank(fits_for(ar1_P5(), y ~ d + x)))

## --- boundary checks of the joint x ---------------------------------------------------------------
bnd_x <- function(d, Nco) { Yc <- matrix(d$y, ncol = max(d$id))[, seq_len(Nco)]; x_ref(Yc, 2) }
boundary <- list(
  wn = data.frame(seed = 1:10, x = vapply(1:10, function(s) bnd_x(ar1_panel(30, 5, 15, 5, rho = 0, seed = 5000 + s), 30), 0)),
  rw30 = data.frame(seed = 1:6, x = vapply(1:6, function(s) bnd_x(ar1_rw_panel(30, 5, 15, 5, 6000 + s), 30), 0)),
  rw100 = data.frame(seed = 1:6, x = vapply(1:6, function(s) bnd_x(ar1_rw_panel(100, 20, 40, 5, 7000 + s), 100), 0)))

## --- R10: the prototype's leave-one-out pick on the motivating design, seeds 1..40 (ar1_mot_panel) ----------
r10 <- data.frame(seed = 1:40, proto_pick = vapply(1:40, function(s) {
  Ym <- matrix(ar1_mot_panel(s)$y, 20, 35)
  unname(pick_row(rank_scores(Ym, 1:30, 31:35, 15, 5))["loo_cell"])
}, 0L))

meta <- list(prototype = "xuyiqing/gsynth-revisit wip/template-v2-2026-10-16 e1eebae", date = "2026-10-16",
             fect_for_fits = as.character(packageVersion("fect")), R = R.version.string)
saveRDS(list(meta = meta, P1 = P1, rules = rules, boundary = boundary, r10 = r10), file.path(gen_dir, "ar1-reference.rds"), version = 2)
cat("wrote", file.path(gen_dir, "ar1-reference.rds"), "\n")
print(P1$path); print(P1$x_joint); print(boundary)
