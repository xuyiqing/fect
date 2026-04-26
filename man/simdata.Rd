\encoding{UTF-8}
\name{simdata}
\docType{data}
\alias{simdata}
\title{Simulated panel data with two latent factors}
\description{
A simulated panel dataset with continuous outcomes used throughout the
package vignettes to demonstrate factor-augmented counterfactual
estimators. The data-generating process follows
\href{https://onlinelibrary.wiley.com/doi/10.1111/ajps.12723}{Liu, Wang,
and Xu (2024)} with one modification (see Format).

The panel has \eqn{N = 200} units and \eqn{T = 35} time periods.
Treatment switches on and off over time (99 of 150 treated units
experience at least one reversal), reflecting a general treatment
pattern rather than simple staggered adoption. The outcome includes two
latent factors (\eqn{r = 2}), so the parallel-trends assumption is
violated and the standard fixed-effects estimator is biased. Treatment
assignment loads on the same factors and fixed effects that enter the
outcome---units with larger \eqn{\lambda_i} and \eqn{\alpha_i} are more
likely to be treated---so the confounding is structural and cannot be
removed by two-way fixed effects alone.
}
\format{
A data frame with the following columns:
\describe{
  \item{id}{unit identifier (1--200)}
  \item{time}{time period (1--35)}
  \item{Y}{observed outcome}
  \item{error}{idiosyncratic error \eqn{\varepsilon_{it} \sim N(0, 2)}}
  \item{eff}{realized treatment effect \eqn{\tau_{it}}}
  \item{tr_cum, tr_prob}{treatment-probability constructions}
  \item{D}{treatment indicator}
  \item{X1, X2}{observed time-varying covariates \eqn{\sim N(0, 1)}
        with coefficients 1 and 3}
  \item{alpha}{unit fixed effect \eqn{\alpha_i \sim N(0, 1)}}
  \item{xi}{time fixed effect \eqn{\xi_t} (AR(1) with drift)}
  \item{F1, F2}{latent time factors \eqn{f_t \in \mathbb{R}^2}
        (one trending, one white noise)}
  \item{L1, L2}{unit-specific factor loadings
        \eqn{\lambda_i \sim N(0.5, 1)}}
  \item{FL1, FL2}{per-cell factor-loading products
        \eqn{\lambda_{i,k} \cdot f_{t,k}} (\eqn{k = 1, 2})}
}

The DGP is
\deqn{Y_{it} = \tau_{it} D_{it} + X_{1,it} + 3 X_{2,it} + \mu
              + 3\alpha_i + \xi_t
              + 2\, \lambda_i' f_t + \varepsilon_{it},}
with grand mean \eqn{\mu = 5} and treatment effect
\eqn{\tau_{it} \sim N(0.4 \cdot \mathrm{tr\_cum}_{it}/T,\; 0.2)}.

The \eqn{2\, \lambda_i' f_t} term doubles the latent factor contribution
relative to the original Liu, Wang, and Xu (2024) DGP. The doubling
strengthens the factor signal-to-noise ratio (variance of the factor
contribution to variance of the residual) from approximately 2.7 to
10.9, which makes the factor structure clearly recoverable by
cross-validated rank-selection procedures on this dataset. The
unmodified DGP is preserved in earlier package versions; see
\code{git log data/simdata.rda} for the prior file.
}
\references{
  Liu, L., Wang, Y., and Xu, Y. (2024).
  A Practical Guide to Counterfactual Estimators for Causal Inference
  with Time-Series Cross-Sectional Data.
  \emph{American Journal of Political Science}, 68(1), 160--176.
}
\keyword{datasets}
