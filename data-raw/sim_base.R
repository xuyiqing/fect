## Generate sim_base: simulated panel data WITHOUT interactive fixed effects
## Same DGP as simdata but with r = 0 (no latent factors)
##
## Model: Y_it = mu + 3*alpha_i + xi_t + X1_it*1 + X2_it*3 + eff_it + error_it
##
## Since treatment assignment in simulateData() depends on factors when r > 0,
## sim_base is constructed from simdata by removing the factor contribution
## from Y and zeroing out F/L/FL columns. This preserves the same units,
## treatment pattern, covariates, and errors as simdata.

load("data/fect.RData")

## Verify consistency: sim_base should equal simdata minus factor structure
stopifnot(all(simdata$D == sim_base$D))
stopifnot(all(simdata$X1 == sim_base$X1))
stopifnot(all(simdata$error == sim_base$error))
stopifnot(all(sim_base$F1 == 0))
stopifnot(all(sim_base$L1 == 0))

usethis::use_data(sim_base, overwrite = TRUE)
