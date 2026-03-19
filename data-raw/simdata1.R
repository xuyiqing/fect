## Generate simdata1: simulated panel data WITHOUT interactive fixed effects
## Same DGP as simdata but with r = 0 (no latent factors)
##
## Model: Y_it = mu + 3*alpha_i + xi_t + X1_it*1 + X2_it*3 + eff_it + error_it
##
## Since treatment assignment in simulateData() depends on factors when r > 0,
## simdata1 is constructed from simdata by removing the factor contribution
## from Y and zeroing out F/L/FL columns. This preserves the same units,
## treatment pattern, covariates, and errors as simdata.

load("data/fect.RData")

## Verify consistency: simdata1 should equal simdata minus factor structure
stopifnot(all(simdata$D == simdata1$D))
stopifnot(all(simdata$X1 == simdata1$X1))
stopifnot(all(simdata$error == simdata1$error))
stopifnot(all(simdata1$F1 == 0))
stopifnot(all(simdata1$L1 == 0))

usethis::use_data(simdata1, overwrite = TRUE)
