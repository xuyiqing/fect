\encoding{UTF-8}
\name{simdata2}
\docType{data}
\alias{simdata2}
\title{Simulated panel data with unit-specific time trends (block DID)}
\description{
A simulated balanced panel with 200 units and 40 time periods.
All 80 treated units receive treatment at period 21 (block DID design).
The DGP includes unit-specific amplitudes on a sinusoidal time trend,
where treated units have systematically larger amplitudes.
}
\format{A data frame with 8000 rows and 4 columns:
\describe{
  \item{id}{Unit identifier (1 to 200)}
  \item{time}{Time period (1 to 40)}
  \item{Y}{Outcome variable}
  \item{D}{Treatment indicator (0/1)}
}
}
\keyword{datasets}
