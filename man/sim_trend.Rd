\encoding{UTF-8}
\name{sim_trend}
\docType{data}
\alias{sim_trend}
\title{Simulated panel data with unit-specific sinusoidal time trends (block DID)}
\description{
A simulated balanced panel with 200 units and 50 time periods.
All 80 treated units receive treatment at period 41 (block DID design).
The DGP includes unit-specific amplitudes on a half-cycle sinusoidal time trend,
where treated units have systematically larger amplitudes.
}
\format{A data frame with 10000 rows and 4 columns:
\describe{
  \item{id}{Unit identifier (1 to 200)}
  \item{time}{Time period (1 to 50)}
  \item{Y}{Outcome variable}
  \item{D}{Treatment indicator (0/1)}
}
}
\keyword{datasets}
