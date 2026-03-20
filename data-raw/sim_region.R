## Generate sim_region: unbalanced panel with region×time interaction fixed effects
## This script produces data/sim_region.rda
##
## DGP: Y_it = alpha_i + xi_t + delta_{g(i),t} + tau * D_it + e_it
## where delta_{g(i),t} are region-specific linear time trends (confounding source).
## Treatment probability and timing depend on region. Units in higher-numbered
## regions enter the panel later, creating unbalancedness.

set.seed(42)
N <- 500; TT <- 20; n.groups <- 5

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

tau <- 1

sim_region <- expand.grid(time = 1:TT, id = 1:N)
sim_region$region <- region[sim_region$id]
sim_region$Y <- alpha_i[sim_region$id] + xi_t[sim_region$time] +
                delta_gt[cbind(region[sim_region$id], sim_region$time)] +
                rnorm(N * TT, 0, 1)
sim_region$D <- 0
for (i in which(is_treated == 1)) {
  idx <- sim_region$id == i & sim_region$time >= treat_time[i]
  sim_region$D[idx] <- 1
}
sim_region$Y <- sim_region$Y + tau * sim_region$D

# Create region × time interaction for group-period fixed effects
sim_region$region_time <- interaction(sim_region$region, sim_region$time, drop = TRUE)

# Staggered panel entry: higher-region units enter later
entry_period <- pmax(1, region[sim_region$id] - 1 +
                     sample(0:1, nrow(sim_region), replace = TRUE))
sim_region <- sim_region[sim_region$time >= entry_period, ]

usethis::use_data(sim_region, overwrite = TRUE)
