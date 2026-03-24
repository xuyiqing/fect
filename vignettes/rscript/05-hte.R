##############################
# 05-hte.R
# Generated from 05-hte.Rmd
##############################
rm(list = ls())
set.seed(1234)

library(fect)
data(sim_base)

## --- hte_setup ---
out.fect <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  method = "fe", force = "two-way", se = TRUE,
  cores = 8, parallel = TRUE, nboots = 200)

## --- hte ---
plot(out.fect, type = "box", xlim = c(-15, 10))

## --- hte_time ---
plot(out.fect, type = "calendar", xlim = c(1, 35))

## --- hte_X1 ---
plot(out.fect, type = "hte", covariate = "X1")

## --- hte_discrete ---
sim_base$X3 <- sample(1:3, size = nrow(sim_base), replace = TRUE)
out.fect.X3 <- fect(Y ~ D + X1 + X2 + X3, data = sim_base, index = c("id","time"),
                   method = "fe", se = TRUE, seed = 123,
                   cores = 8, nboots = 200, parallel = TRUE)

## --- plot-hte-discrete ---
plot(out.fect.X3, type="hte", covariate = "X3",
     xlab = "", ylab = "Effect of D on Y",
     covariate.labels = c("USA", "China", "UK"),
     ylim = c(-2, 6))

## --- cm_fe ---
out.cm <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id", "time"),
               method = "fe", force = "two-way", se = TRUE,
               cm = TRUE, parallel = TRUE, cores = 8, nboots = 200)

## --- hte_em ---
plot(out.cm, type = "hte", covariate = "X1",
     xlab = "Moderator (X1)", ylab = "Effect on Y")

## --- hte_cm ---
plot(out.cm, type = "hte", covariate = "X1", cm = TRUE,
     xlab = "Moderator (X1)", ylab = "Effect on Y")

## --- hte_scatter ---
plot(out.cm, type = "hte", covariate = "X1", cm = TRUE, loess.fit = FALSE,
     xlab = "Moderator (X1)", ylab = "Effect on Y")

## --- hte_placebo ---
plot(out.cm, type = "hte", covariate = "X1",
     pretreatment = TRUE, num.pretreatment = 3,
     xlab = "X1", ylab = "Placebo Effect")

## --- hte_placebo_cm ---
plot(out.cm, type = "hte", covariate = "X1", cm = TRUE,
     pretreatment = TRUE, num.pretreatment = 3,
     xlab = "X1", ylab = "Placebo Effect")

## --- iden_test ---
iden.test <- fect_iden(out.cm, moderator = "X1")

## --- iden_results ---
cat("=== Treated cells (e1) ===\n")
cat("  n =", iden.test$e1$n, "\n")
cat("  R-squared =", round(iden.test$e1$r2, 4), "\n")
cat("  Test stat =", round(iden.test$e1$stat, 3), "\n")
cat("  df =", iden.test$e1$df, "\n")
cat("  p-value =", round(iden.test$e1$p, 4), "\n\n")
cat("=== Control cells (e0) ===\n")
cat("  n =", iden.test$e0$n, "\n")
cat("  R-squared =", round(iden.test$e0$r2, 4), "\n")
cat("  Test stat =", round(iden.test$e0$stat, 3), "\n")
cat("  df =", iden.test$e0$df, "\n")
cat("  p-value =", round(iden.test$e0$p, 4), "\n")

## --- iden_components ---
iden.quad <- fect_iden(out.cm, moderator = "X1", interaction = FALSE)
cat("Quadratic-only: p =",
    round(iden.quad$e1$p, 4), "(treated),",
    round(iden.quad$e0$p, 4), "(control)\n")
iden.inter <- fect_iden(out.cm, moderator = "X1", quadratic = FALSE)
cat("Interaction-only: p =",
    round(iden.inter$e1$p, 4), "(treated),",
    round(iden.inter$e0$p, 4), "(control)\n")

## --- discrete_cm_setup ---
sim_base$X3 <- sample(1:3, size = nrow(sim_base), replace = TRUE)
out.discrete <- fect(Y ~ D + X1 + X2 + X3, data = sim_base,
                     index = c("id", "time"),
                     method = "fe", force = "two-way", se = TRUE,
                     cm = TRUE, parallel = TRUE, cores = 8, nboots = 200)

## --- discrete_em_plot ---
plot(out.discrete, type = "hte", covariate = "X3",
     covariate.labels = c("USA", "China", "UK"),
     xlab = "", ylab = "Effect on Y")

## --- discrete_cm_plot ---
plot(out.discrete, type = "hte", covariate = "X3", cm = TRUE,
     covariate.labels = c("USA", "China", "UK"),
     xlab = "", ylab = "Effect on Y")
