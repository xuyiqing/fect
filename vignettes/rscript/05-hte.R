##############################
# 05-hte.R
# Generated from 05-hte.Rmd
##############################
rm(list = ls())
set.seed(1234)

library(fect)
data(fect)

## --- hte_setup ---
out.fect <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  method = "fe", force = "two-way", se = TRUE,
  cores = 8, parallel = TRUE, nboots = 1000)

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
                   cores = 8, nboots = 1000, parallel = TRUE)

## --- plot-hte-discrete ---
plot(out.fect.X3, type="hte", covariate = "X3",
     xlab = "", ylab = "Effet of D on Y",
     covariate.labels = c("USA", "China", "UK"),
     ylim = c(-2, 6))
