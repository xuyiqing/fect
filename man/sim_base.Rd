\encoding{UTF-8}
\name{sim_base}
\docType{data}
\alias{sim_base}
\title{Simulated data (no interactive fixed effects)}
\description{A simulated panel with continuous outcomes and no latent factors, so that the two-way fixed effects model is correctly specified. It shares the treatment, the covariates, the unit and time effects, the errors and the effects with \code{\link{simdata}}; only the outcome differs. Used in the chapter on the imputation estimator of the user manual.

The panel is balanced, with 200 units (IDs 101 to 300) and 35 periods. 150 units are treated at some point and 50 never; 99 treated units have at least one treatment reversal.}
\format{A data frame with 7,000 rows and 18 columns:
\describe{
  \item{id}{unit identifier (101 to 300).}
  \item{time}{time period (1 to 35).}
  \item{Y}{observed outcome, \eqn{Y_{it} = \tau_{it} D_{it} + X_{1,it} + 3 X_{2,it} + \alpha_i + \xi_t + \varepsilon_{it}}{Y = eff * D + X1 + 3 * X2 + alpha + xi + error}. Unlike \code{simdata}, the outcome has no factor term and no grand mean, and the unit effect enters with coefficient 1.}
  \item{error}{idiosyncratic error \eqn{\varepsilon_{it}}, as in \code{simdata}.}
  \item{eff}{realized treatment effect \eqn{\tau_{it}}; 0 in untreated cells.}
  \item{tr_cum, tr_prob}{as in \code{simdata}: \code{tr_cum} counts the consecutive treated periods up to and including the current one (0 when untreated), and \code{tr_prob} is the treatment probability (between 0 and 1).}
  \item{D}{treatment indicator (0/1).}
  \item{X1, X2}{observed time-varying covariates, with coefficients 1 and 3.}
  \item{alpha}{unit effect \eqn{\alpha_i}.}
  \item{xi}{time effect \eqn{\xi_t}.}
  \item{F1, F2, L1, L2, FL1, FL2}{the factor columns of \code{simdata}, set to 0 here.}
}
}
\seealso{\code{\link{simdata}}}
\keyword{datasets}
\references{
  Liu, L., Wang, Y., and Xu, Y. (2024).
  A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional Data.
  \emph{American Journal of Political Science, 68}(1), 160-176.
}
