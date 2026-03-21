##########################
# Install Packages
##########################
rm(list = ls())
# install.packages("fect")
# devtools::install_github("xuyiqing/fect")
library(fect)
data(fect)

##########################
# From: ./05-hte.Rmd
##########################

set.seed(1234)
rm(list = ls())
library(fect)
data(fect)

##########################
# Causal Moderation (FE)
##########################

out.cm <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
               method = "fe", force = "two-way", se = TRUE,
               cm = TRUE, parallel = TRUE, cores = 4, nboots = 200)

## Standard HTE plot (descriptive)
plot(out.cm, type = "hte", covariate = "X1")

## Causal moderation plot (model-based)
plot(out.cm, type = "hte", covariate = "X1", cm = TRUE)

## Without loess smoothing
plot(out.cm, type = "hte", covariate = "X1", cm = TRUE, loess.fit = FALSE)

##########################
# Over-Identification Test
##########################

iden.test <- fect_iden(out.cm, moderator = "X1")

## Treated cells (e1)
cat("Treated cells test:\n")
cat("  n =", iden.test$e1$n, "\n")
cat("  R^2 =", round(iden.test$e1$r2, 4), "\n")
cat("  Test stat =", round(iden.test$e1$stat, 3), "\n")
cat("  df =", iden.test$e1$df, "\n")
cat("  p-value =", round(iden.test$e1$p, 4), "\n")

## Control cells (e0)
cat("\nControl cells test:\n")
cat("  n =", iden.test$e0$n, "\n")
cat("  R^2 =", round(iden.test$e0$r2, 4), "\n")
cat("  Test stat =", round(iden.test$e0$stat, 3), "\n")
cat("  df =", iden.test$e0$df, "\n")
cat("  p-value =", round(iden.test$e0$p, 4), "\n")

##########################
# Causal Moderation (IFE)
##########################

out.cm.ife <- fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                   method = "ife", force = "two-way", CV = TRUE, r = c(0, 5),
                   cm = TRUE, se = TRUE, parallel = TRUE, cores = 4, nboots = 200)

plot(out.cm.ife, type = "hte", covariate = "X1", cm = TRUE)

iden.test.ife <- fect_iden(out.cm.ife, moderator = "X1")
cat("IFE - Treated cells p-value:", round(iden.test.ife$e1$p, 4), "\n")
cat("IFE - Control cells p-value:", round(iden.test.ife$e0$p, 4), "\n")
