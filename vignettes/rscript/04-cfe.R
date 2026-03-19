##########################
# Complex Fixed Effects (CFE)
##########################
rm(list = ls())
library(fect)
data(fect)

set.seed(1234)

##########################
# 4.2 Additional fixed effects
##########################

## DGP
set.seed(42)
N <- 300; TT <- 20; n.groups <- 5

region <- rep(1:n.groups, each = N / n.groups)
alpha_i <- rnorm(N, 0, 1)
xi_t <- rnorm(TT, 0, 0.5)

region_slope <- c(-1, -0.5, 0, 0.5, 1) * 0.08
delta_gt <- outer(region_slope, (1:TT) - TT / 2)

treat_prob <- c(0.05, 0.15, 0.30, 0.60, 0.80)
is_treated <- rbinom(N, 1, treat_prob[region])
treat_time <- rep(NA, N)
for (i in which(is_treated == 1)) {
  treat_time[i] <- max(5, min(16,
    round(13 - 2 * (region[i] - 3)) + sample(-1:1, 1)))
}

tau <- 2

dat_fe <- expand.grid(time = 1:TT, id = 1:N)
dat_fe$region <- region[dat_fe$id]
dat_fe$Y <- alpha_i[dat_fe$id] + xi_t[dat_fe$time] +
            delta_gt[cbind(region[dat_fe$id], dat_fe$time)] +
            rnorm(N * TT, 0, 1)
dat_fe$D <- 0
for (i in which(is_treated == 1)) {
  idx <- dat_fe$id == i & dat_fe$time >= treat_time[i]
  dat_fe$D[idx] <- 1
}
dat_fe$Y <- dat_fe$Y + tau * dat_fe$D

# Staggered panel entry: higher-region units enter later
entry_period <- pmax(1, region[dat_fe$id] - 1 +
                     sample(0:1, nrow(dat_fe), replace = TRUE))
dat_fe <- dat_fe[dat_fe$time >= entry_period, ]

## FE only
out.fe.only <- fect(Y ~ D, data = dat_fe,
  index = c("id", "time"),
  method = "fe", force = "two-way",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.fe.only, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "FE Only — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## With region FE
out.cfe.region <- fect(Y ~ D, data = dat_fe,
  index = c("id", "time", "region"),
  method = "cfe", force = "two-way",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.cfe.region, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with Region FE — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

##########################
# 4.3 Time-invariant covariates with time-varying coefficients
##########################

## FE baseline (simdata)
out.fe.base <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "fe", force = "two-way",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.fe.base, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "FE Only (simdata) — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## With Z = L1
simdata$gamma_t <- simdata$time

out.cfe.z <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "L1", gamma = "gamma_t",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.cfe.z, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with Z = L1 — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

##########################
# 4.4 Unit-specific time trends
##########################

## DGP
set.seed(42)
N2 <- 200; TT2 <- 40; T0 <- 21

alpha2_i <- rnorm(N2, 0, 1)
xi2_t <- rnorm(TT2, 0, 0.5)
is_treated2 <- c(rep(1, N2 * 0.4), rep(0, N2 * 0.6))
kappa2_i <- ifelse(is_treated2 == 1,
                   runif(N2, 2, 4),
                   runif(N2, 0.5, 1.5))
Q2_t <- sin(2 * pi * (1:TT2) / TT2)
tau2 <- 3

simdata2 <- expand.grid(time = 1:TT2, id = 1:N2)
simdata2$Y <- alpha2_i[simdata2$id] + xi2_t[simdata2$time] +
              kappa2_i[simdata2$id] * Q2_t[simdata2$time] +
              rnorm(N2 * TT2, 0, 1)
simdata2$D <- as.integer(is_treated2[simdata2$id] == 1 & simdata2$time >= T0)
simdata2$Y <- simdata2$Y + tau2 * simdata2$D

## panelView
library(panelView)
panelview(Y ~ D, data = simdata2, index = c("id", "time"),
  axis.lab = "time", xlab = "Time", ylab = "Unit",
  gridOff = TRUE, by.timing = TRUE,
  background = "white", main = "simdata2: Treatment Status")

## FE only
out.fe.trend <- fect(Y ~ D, data = simdata2,
  index = c("id", "time"),
  method = "fe", force = "two-way",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.fe.trend, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "FE Only (simdata2) — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## Linear trend
out.cfe.lin <- fect(Y ~ D, data = simdata2,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Q.type = "linear",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.cfe.lin, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with Linear Trend — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

## B-spline trend
out.cfe.bs <- fect(Y ~ D, data = simdata2,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Q.type = "bspline",
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.cfe.bs, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE with B-spline Trend — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)

##########################
# 4.5 CFE with factors
##########################
simdata$gamma_t <- simdata$time

## Model comparison
out.fe <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "fe", force = "two-way", se = FALSE)

out.cfe.z.only <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "L1", gamma = "gamma_t",
  se = FALSE)

out.cfe.z.f1 <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "L1", gamma = "gamma_t",
  r = 1, se = FALSE)

out.ife.r2 <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "ife", force = "two-way",
  r = 2, se = FALSE)

mspe.out <- fect_mspe(
  list(FE = out.fe,
       CFE_Z = out.cfe.z.only,
       CFE_Z_F1 = out.cfe.z.f1,
       IFE_r2 = out.ife.r2),
  seed = 1234)
print(mspe.out$summary[, c("Model", "MSPE")])

## Best model placebo test
out.cfe.best <- fect(Y ~ D + X1 + X2, data = simdata,
  index = c("id", "time"),
  method = "cfe", force = "two-way",
  Z = "L1", gamma = "gamma_t",
  r = 1,
  se = TRUE, parallel = TRUE, nboots = 200,
  placeboTest = TRUE, placebo.period = c(-2, 0))
plot(out.cfe.best, cex.text = 0.8,
     stats = c("placebo.p", "equiv.p"),
     main = "CFE (Z + 1 Factor) — Placebo Test",
     cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8)
