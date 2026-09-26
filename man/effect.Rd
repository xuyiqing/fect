\name{effect}
\alias{effect}
\title{Calculate Cumulative or Sub-group Treatment Effects}
\description{
  Calculates cumulative or average treatment effects for specified units and time periods based on a fitted \code{fect} object. The function supports both cumulative effects over time and period-specific average treatment effects, with bootstrap-based uncertainty estimates.
}

\usage{
effect(
  x,
  cumu = TRUE,
  id = NULL,
  period = NULL,
  plot = FALSE,
  count = TRUE,
  xlab = NULL,
  ylab = NULL,
  main = NULL
)
}

\arguments{
  \item{x}{A \code{fect} object containing treatment effect estimates and bootstrap results.}

  \item{cumu}{Logical. If \code{TRUE} (default), calculates cumulative treatment effects: the running sum of the per-period average effects. If \code{FALSE}, calculates period-specific average treatment effects.}

  \item{id}{Character vector or NULL. Unit identifiers to include in the analysis; one unit is fine. If \code{NULL} (default), all treated units are included.}

  \item{period}{Numeric vector of length 2 specifying the time window \code{c(start, end)} for effect calculation. If \code{NULL}, uses the maximum possible window based on the data.}

  \item{plot}{Logical. If \code{TRUE}, creates a visualization of the cumulative treatment effects with confidence intervals and a bar chart showing the number of treated units at each time point. Default is \code{FALSE}.}

  \item{count}{Logical. If \code{TRUE}, shows the count bars in the plot.}

  \item{xlab}{Character. X-axis label for the plot.}

  \item{ylab}{Character. Y-axis label for the plot.}

  \item{main}{Character. Main title for the plot.}
}

\details{
  The function processes treatment effects in several steps:

  1. Selects units based on the \code{id} parameter or includes all treated units if \code{id = NULL}.

  2. Calculates relative time to treatment for each unit.

  3. If \code{cumu = TRUE}, computes cumulative effects by summing average effects up to each period. Each period's average effect is the mean over the selected treated units' cells at that event time, so the cumulative effect at event time \eqn{k} is the sum of these averages from \code{period[1]} to \eqn{k}. \code{att.cumu(x, weighted = TRUE)} gives the count-weighted version that \code{fect} used before version 2.4.7; it differs when the number of treated cells changes across event times.

  4. Performs bootstrap analysis to estimate uncertainty (standard errors, confidence intervals, and p-values).

  The function supports different inference methods (bootstrap, jackknife, parametric) and adjusts calculations accordingly. The method is read from the fitted object (\code{x$vartype}), not from the call. For bootstrap fits, the confidence interval is formed by the 2.5\% and 97.5\% quantiles of the replicate estimates and the p-value from the share of replicates on each side of zero. For parametric and jackknife fits (the jackknife SE is scaled by \eqn{\sqrt{N-1}}{sqrt(N - 1)}), the interval is the estimate plus or minus a t critical value times the SE, with a t-based p-value; parametric replicates are simulated under no treatment effect, so their quantiles are not an interval for the estimate.

  Note: The function requires bootstrap results in the input \code{fect} object (\code{keep.sims = TRUE} must be set when fitting the model). It returns \code{NULL}, with a warning, for panels with treatment reversals.
}

\value{
  Returns the input \code{fect} object with two added elements:
  \item{effect.est.avg}{Vector of point estimates for cumulative or average treatment effects, one per period in \code{period}.}
  \item{effect.est.att}{Matrix with one row per period and the following columns:
    \itemize{
      \item ATT: Point estimates
      \item S.E.: Standard errors
      \item CI.lower: Lower bound of confidence interval
      \item CI.upper: Upper bound of confidence interval
      \item p.value: Two-sided p-values
    }
  }
}

\section{Warning}{
  The function will stop with an error if:
  \itemize{
    \item No bootstrap results are available in the input object
    \item The specified ending period exceeds the maximum available period
  }
  For a panel with treatment reversals it warns and returns \code{NULL}.
}

\examples{
\donttest{
library(fect)
data(simgsynth) # no treatment reversals
fit <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
            method = "fe", force = "two-way", se = TRUE, nboots = 50,
            parallel = FALSE, keep.sims = TRUE)

# Cumulative effects for all treated units
results <- effect(fit)
head(results$effect.est.att)

# Period-specific effects for two treated units over event times 1 to 4
ids <- head(unique(simgsynth$id[simgsynth$D == 1]), 2)
results_specific <- effect(fit, cumu = FALSE, id = ids, period = c(1, 4))
results_specific$effect.est.att
}
}

\seealso{
  \code{\link{fect}}, \code{\link{plot.fect}}
}

\references{
  Liu, L., Wang, Y., & Xu, Y. (2022).
  A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional Data.
  \emph{American Journal of Political Science}, 68(1), 160-176.
}

\author{
  Shiyun Hu, Licheng Liu, Ye Wang, and Yiqing Xu
}
