# Full conformal vs Jackknife+ for target-date rules (2026-11-11)

Script: `tests/coverage-study/conformal_full_coverage.R` (2,000 replications, 8 cores). One treated unit, one
target date, T0 = 20 pre-treatment dates, effect 2, gamma = 0.1. Designs: `additive` (unit and time effects plus
iid noise: units exchangeable), `factor_in` (two factors, loadings uniform on [0, 1] for every unit),
`factor_out` (the treated unit's loadings at 1.6, outside the controls' range). Rules: `did` and ridge-simplex
`sc` (eta = 1e-5). Guaranteed coverage: full conformal 1 - floor(0.1 N)/N = 0.90 at N = 20 and N = 40;
Jackknife+ 1 - (2q - 1)/N = 0.85 at N = 20 and 0.875 at N = 40.

|  J | design     | rule | full conformal | Jackknife+ | FC width | J+ width |
|---:|------------|------|---------------:|-----------:|---------:|---------:|
| 19 | additive   | did  | 0.901 | 0.941 | 3.60 | 4.39 |
| 19 | additive   | sc   | 0.902 | 0.949 | 4.08 | 5.04 |
| 19 | factor_in  | did  | 0.897 | 0.949 | 3.64 | 4.35 |
| 19 | factor_in  | sc   | 0.893 | 0.955 | 2.36 | 2.94 |
| 19 | factor_out | did  | 0.334 | 0.376 | 3.66 | 4.38 |
| 19 | factor_out | sc   | 0.246 | 0.325 | 2.60 | 3.11 |
| 39 | additive   | did  | 0.881 | 0.948 | 3.44 | 4.19 |
| 39 | additive   | sc   | 0.887 | 0.944 | 3.84 | 4.66 |
| 39 | factor_in  | did  | 0.904 | 0.955 | 3.47 | 4.10 |
| 39 | factor_in  | sc   | 0.905 | 0.955 | 2.12 | 2.57 |
| 39 | factor_out | did  | 0.318 | 0.359 | 3.43 | 4.07 |
| 39 | factor_out | sc   | 0.241 | 0.297 | 2.31 | 2.60 |

Monte Carlo s.e. of a coverage near 0.9 is 0.0067. Full conformal sets sit at their guarantee and are 15-20%
narrower than Jackknife+ intervals, which over-cover (designs where both hold). No set was disconnected. When the treated unit lies outside
the controls' range, both fail, as the theory allows (the scores are not exchangeable).

The J = 39 additive DID cell (0.881) is 2.8 s.e. below 0.90. A rerun with 20,000 replications gave 0.896; the same
seeds in plain base R (ranking |D_i - mean(D)|, which is what the DID scores reduce to) give the same 0.104 miss
rate, and 200,000 fresh iid draws give 0.1003. So the shortfall is Monte Carlo noise of those seeds, not the code.
