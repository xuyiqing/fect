##########################
# Effect Heterogeneity
##########################
rm(list = ls())
library(fect)
data(fect)

set.seed(1234)

##########################
# From: ./05-hte.Rmd
##########################

##########################
# Setup
##########################

out.fect <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id","time"),
  method = "fe", force = "two-way", se = TRUE,
  cores = 8, parallel = TRUE, nboots = 1000)

##########################
# Box plot
##########################

plot(out.fect, type = "box", xlim = c(-15, 10))

##########################
# By calendar time
##########################

plot(out.fect, type = "calendar", xlim = c(1, 35))

##########################
# By a covariate
##########################

plot(out.fect, type = "hte", covariate = "X1")

# Discrete covariate
sim_base$X3 <- sample(1:3, size = nrow(sim_base), replace = TRUE)
out.fect.X3 <- fect(Y ~ D + X1 + X2 + X3, data = sim_base, index = c("id","time"),
                   method = "fe", se = TRUE, seed = 123,
                   cores = 8, nboots = 1000, parallel = TRUE)

plot(out.fect.X3, type="hte", covariate = "X3",
     xlab = "", ylab = "Effet of D on Y",
     covariate.labels = c("USA", "China", "UK"),
     ylim = c(-2, 6))
