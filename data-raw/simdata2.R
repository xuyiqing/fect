## Generate simdata2: block DID with unit-specific sinusoidal time trends
## This script produces data/simdata2.rda

set.seed(42)
N2 <- 200; TT2 <- 40; T0 <- 21  # treatment starts at period 21

# Unit-level
alpha2_i <- rnorm(N2, 0, 1)

# Time-level
xi2_t <- rnorm(TT2, 0, 0.5)

# Unit-specific trend amplitude
# kappa_i ~ Uniform(0.5, 1.5) for control, Uniform(2, 4) for treated
is_treated2 <- c(rep(1, N2 * 0.4), rep(0, N2 * 0.6))  # 80 treated, 120 control
kappa2_i <- ifelse(is_treated2 == 1,
                   runif(N2, 2, 4),
                   runif(N2, 0.5, 1.5))

# True time trend: sin wave
Q2_t <- sin(2 * pi * (1:TT2) / TT2)

# Treatment effect
tau2 <- 3

# Build panel
simdata2 <- expand.grid(time = 1:TT2, id = 1:N2)
simdata2$Y <- alpha2_i[simdata2$id] + xi2_t[simdata2$time] +
              kappa2_i[simdata2$id] * Q2_t[simdata2$time] +
              rnorm(N2 * TT2, 0, 1)
simdata2$D <- as.integer(is_treated2[simdata2$id] == 1 & simdata2$time >= T0)
simdata2$Y <- simdata2$Y + tau2 * simdata2$D

usethis::use_data(simdata2, overwrite = TRUE)
