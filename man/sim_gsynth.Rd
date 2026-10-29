\encoding{UTF-8}
\name{sim_gsynth}
\docType{data}
\alias{sim_gsynth}
\title{Simulated data for Gsynth (no treatment reversal)}
\description{A simulated panel with no treatment reversal, used to demonstrate the generalized synthetic control method (the never-treated estimation regime) in the user manual.

The panel is balanced, with 50 units (IDs 101 to 150) and 30 periods. Units 101 to 105 are treated from period 21 to the end; the other 45 units are never treated. The outcome has two latent factors.}
\format{A data frame with 1,500 rows and 15 columns:
\describe{
  \item{id}{unit identifier (101 to 150).}
  \item{time}{time period (1 to 30).}
  \item{Y}{observed outcome, \eqn{Y_{it} = \tau_{it} D_{it} + X_{1,it} + 3 X_{2,it} + \mu + \alpha_i + \xi_t + F_{1t} L_{1i} + F_{2t} L_{2i} + \varepsilon_{it}}{Y = eff * D + X1 + 3 * X2 + mu + alpha + xi + F1 * L1 + F2 * L2 + error}.}
  \item{D}{treatment indicator (0/1).}
  \item{X1, X2}{observed time-varying covariates, with coefficients 1 and 3.}
  \item{eff}{treatment effect \eqn{\tau_{it}}, for every unit: 0 up to period 20, and from period 21 an effect whose average over units rises from about 0.7 to about 9.7 by period 30. It enters \code{Y} only where \code{D = 1}.}
  \item{error}{idiosyncratic error \eqn{\varepsilon_{it}}.}
  \item{mu}{grand mean (5).}
  \item{alpha}{unit effect \eqn{\alpha_i}.}
  \item{xi}{time effect \eqn{\xi_t}.}
  \item{F1, F2}{latent time factors.}
  \item{L1, L2}{unit-specific factor loadings.}
}
}
\seealso{\code{\link{simgsynth}}, an identical copy.}
\keyword{datasets}
\references{
  Xu, Y. (2017).
  Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models.
  \emph{Political Analysis, 25}(1), 57-76.
}
