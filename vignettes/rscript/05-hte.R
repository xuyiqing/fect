## ----.common, include = FALSE-------------------------------------------------
source("_common.R")


## ----setup-hte, echo = FALSE, message = FALSE, warning = FALSE----------------
set.seed(1234)
data(sim_base)


## ----hte_setup, eval=TRUE, cache=TRUE, message=FALSE, results='hide'----------
out.fect <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  method = "fe", force = "two-way", se = TRUE,
  parallel = TRUE, cores = 16, nboots = 1000)


## ----hte, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.fect, type = "box", xlim = c(-15, 10))


## ----hte_time, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.fect, type = "calendar", xlim = c(1, 35))


## ----hte_X1, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.fect, type = "hte", covariate = "X1")


## ----hte_discrete, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
sim_base$X3 <- sample(1:3, size = nrow(sim_base), replace = TRUE)
out.fect.X3 <- fect(Y ~ D + X1 + X2 + X3, data = sim_base, index = c("id","time"),
                   method = "fe", se = TRUE, seed = 123,
                   nboots = 1000, parallel = TRUE, cores = 16)


## ----plot-hte-discrete, fig.width = 6, fig.height = 4.5-----------------------
plot(out.fect.X3, type="hte", covariate = "X3",
     xlab = "", ylab = "Effect of D on Y",
     covariate.labels = c("USA", "China", "UK"),
     ylim = c(-2, 6))


## ----cm_fe, eval = TRUE, cache = TRUE, message = FALSE, results = 'hide'------
out.cm <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id", "time"),
               method = "fe", force = "two-way", se = TRUE,
               cm = TRUE, parallel = TRUE, cores = 16, nboots = 1000)


## ----hte_em, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.cm, type = "hte", covariate = "X1",
     xlab = "Moderator (X1)", ylab = "Effect on Y")


## ----hte_cm, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.cm, type = "hte", covariate = "X1", cm = TRUE,
     xlab = "Moderator (X1)", ylab = "Effect on Y", ylim = c(-0.5, 5))


## ----hte_scatter, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.cm, type = "hte", covariate = "X1", cm = TRUE, loess.fit = FALSE,
     xlab = "Moderator (X1)", ylab = "Effect on Y", ylim = c(-0.5, 5))


## ----hte_placebo, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.cm, type = "hte", covariate = "X1",
     pretreatment = TRUE, num.pretreatment = 3,
     xlab = "X1", ylab = "Placebo Effect")


## ----hte_placebo_cm, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.cm, type = "hte", covariate = "X1", cm = TRUE,
     pretreatment = TRUE, num.pretreatment = 3,
     xlab = "X1", ylab = "Placebo Effect", ylim = c(-0.5, 1.5))


## ----iden_test, eval = TRUE, cache = TRUE, message = FALSE--------------------
iden.test <- fect_iden(out.cm, moderator = "X1")


## ----iden_results, eval = TRUE------------------------------------------------
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


## ----iden_components, eval = TRUE, cache = TRUE, message = FALSE--------------
# Quadratic terms only (no interactions)
iden.quad <- fect_iden(out.cm, moderator = "X1", interaction = FALSE)
cat("Quadratic-only: p =",
    round(iden.quad$e1$p, 4), "(treated),",
    round(iden.quad$e0$p, 4), "(control)\n")

# Interactions only (no quadratics)
iden.inter <- fect_iden(out.cm, moderator = "X1", quadratic = FALSE)
cat("Interaction-only: p =",
    round(iden.inter$e1$p, 4), "(treated),",
    round(iden.inter$e0$p, 4), "(control)\n")


## ----discrete_cm_setup, eval = TRUE, cache = TRUE, message = FALSE, results = 'hide'----
sim_base$X3 <- sample(1:3, size = nrow(sim_base), replace = TRUE)
out.discrete <- fect(Y ~ D + X1 + X2 + X3, data = sim_base,
                     index = c("id", "time"),
                     method = "fe", force = "two-way", se = TRUE,
                     cm = TRUE, parallel = TRUE, cores = 16, nboots = 1000)


## ----discrete_em_plot, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.discrete, type = "hte", covariate = "X3",
     covariate.labels = c("USA", "China", "UK"),
     xlab = "", ylab = "Effect on Y", ylim = c(-0.5, 5))


## ----discrete_cm_plot, eval = TRUE, cache = TRUE, warning = FALSE, fig.width = 6, fig.height = 4.5----
plot(out.discrete, type = "hte", covariate = "X3", cm = TRUE,
     covariate.labels = c("USA", "China", "UK"),
     xlab = "", ylab = "Effect on Y", ylim = c(-0.5, 5))

