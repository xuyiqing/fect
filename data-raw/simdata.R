## Generate simdata: simulated panel data with interactive fixed effects
## DGP from Liu, Wang, and Xu (2022, AJPS)
## Original source: LWX replication code, 1_ex_sim0.R + simulateData.R
##
## Model: Y_it = mu + 3*alpha_i + xi_t + X1_it*1 + X2_it*3 + L_i'F_t + eff_it + error_it
##
## Key features:
##   N = 200 (150 ever-treated, 50 never-treated)
##   T = 35
##   Staggered adoption with treatment reversals
##   2 latent factors: F1 (trend type), F2 (white noise)
##   Factor sizes: Fsize = c(1.5, 1)
##   Treatment effect: eff ~ N(0.4 * tr_cum/T, 0.2) * D
##   Error sd = 2, treatment noise = 0.5
##   mu = 5 (grand mean)

source("data-raw/simulateData.R")

set.seed(12345)
simdata <- simulateData(
  N = 200, TT = 35, r = 2,
  p = 2, beta = c(1, 3),
  force = 3, mu = 5,
  Rtype = "n",
  Ftype = c("trend", "white"),
  Fsize = c(1.5, 1),
  eff.size = 0.4,
  eff.noise = 0.2,
  seed = 12345
)

usethis::use_data(simdata, overwrite = TRUE)
