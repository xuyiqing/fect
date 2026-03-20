## Generate sim_trend: block DID with unit-specific sinusoidal time trends
## This script produces data/sim_trend.rda
##
## DGP: Y_it = alpha_i + xi_t + kappa_i * sin(2*pi*t / (2*TT)) + tau * D_it + e_it
## where kappa_i ~ U(0.5, 1.0) for treated, U(0.125, 0.375) for controls.
## N = 200 (80 treated, 120 control), T = 50, T0 = 41.

set.seed(42)
N2 <- 200; TT2 <- 50; T0 <- 41

alpha2_i <- rnorm(N2, 0, 1)
xi2_t <- rnorm(TT2, 0, 0.5)

is_treated2 <- c(rep(1, N2 * 0.4), rep(0, N2 * 0.6))
kappa2_i <- ifelse(is_treated2 == 1,
                   runif(N2, 0.5, 1.0),
                   runif(N2, 0.125, 0.375))

# True time trend: half-cycle sin wave
Q2_t <- sin(2 * pi * (1:TT2) / (2 * TT2))

tau2 <- 1

sim_trend <- expand.grid(time = 1:TT2, id = 1:N2)
sim_trend$Y <- alpha2_i[sim_trend$id] + xi2_t[sim_trend$time] +
               kappa2_i[sim_trend$id] * Q2_t[sim_trend$time] +
               rnorm(N2 * TT2, 0, 1)
sim_trend$D <- as.integer(is_treated2[sim_trend$id] == 1 & sim_trend$time >= T0)
sim_trend$Y <- sim_trend$Y + tau2 * sim_trend$D

usethis::use_data(sim_trend, overwrite = TRUE)
