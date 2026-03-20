## Generate sim_linear: block DID with unit-specific linear time trends
## This script produces data/sim_linear.rda
##
## DGP: Y_it = alpha_i + xi_t + kappa_i * (t/T) + tau * D_it + e_it
## where kappa_i ~ U(0.5, 1.0) for treated, U(0.1, 0.3) for controls.
## N = 200 (80 treated, 120 control), T = 50, T0 = 41.

set.seed(42)
N_lin <- 200; TT_lin <- 50; T0_lin <- 41

alpha_lin_i <- rnorm(N_lin, 0, 1)
xi_lin_t <- rnorm(TT_lin, 0, 0.5)

is_treated_lin <- c(rep(1, N_lin * 0.4), rep(0, N_lin * 0.6))
kappa_lin_i <- ifelse(is_treated_lin == 1,
                      runif(N_lin, 0.5, 1.0),
                      runif(N_lin, 0.1, 0.3))

# True trend: linear t/TT with unit-specific slope
Q_lin_t <- (1:TT_lin) / TT_lin

tau_lin <- 1

sim_linear <- expand.grid(time = 1:TT_lin, id = 1:N_lin)
sim_linear$Y <- alpha_lin_i[sim_linear$id] + xi_lin_t[sim_linear$time] +
                kappa_lin_i[sim_linear$id] * Q_lin_t[sim_linear$time] +
                rnorm(N_lin * TT_lin, 0, 1)
sim_linear$D <- as.integer(is_treated_lin[sim_linear$id] == 1 &
                            sim_linear$time >= T0_lin)
sim_linear$Y <- sim_linear$Y + tau_lin * sim_linear$D

usethis::use_data(sim_linear, overwrite = TRUE)
