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

  \item{cumu}{Logical. If \code{TRUE} (default), calculates cumulative treatment effects. If \code{FALSE}, calculates period-specific average treatment effects.}

  \item{id}{Character vector or NULL. Unit identifiers to include in the analysis. If \code{NULL} (default), all treated units are included.}

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

  3. If \code{cumu = TRUE}, computes cumulative effects by summing average effects up to each period.

  4. Performs bootstrap analysis to estimate uncertainty (standard errors, confidence intervals, and p-values).

  The function supports different inference methods (bootstrap, jackknife, parametric) and adjusts calculations accordingly.

  Note: The function requires bootstrap results in the input \code{fect} object (\code{keep.sims = TRUE} must be set when fitting the model).
}

\value{
  Returns a list containing:
  \item{eff}{Vector of point estimates for cumulative or average treatment effects.}
  \item{est.eff}{Matrix containing the following columns:
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
    \item The panel contains treatment reversals
    \item The specified ending period exceeds the maximum available period
  }
}

\examples{
\dontrun{
# Fit fect model with bootstrap
fit <- fect(Y ~ D + X, data = panel_data, keep.sims = TRUE)

# Calculate cumulative effects for all treated units
results <- effect(fit)

# Calculate period-specific effects for specific units
results_specific <- effect(fit,
                          cumu = FALSE,
                          id = c("unit1", "unit2"),
                          period = c(1, 4))

# View results
print(results$est.catt)
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
