##########################
# Install Packages 
##########################
rm(list = ls())
install.packages("fect")
devtools::install_github("xuyiqing/fect")
installed.packages()["fect", "Version"]
devtools::install_github('xuyiqing/panelView')
install_all <- function(packages) {
  installed_pkgs <- installed.packages()[, "Package"]
  for (pkg in packages) {
    if (!pkg %in% installed_pkgs) {
      install.packages(pkg)
    }
  }
}
packages <- c("abind", "doParallel", "doRNG", "fixest", "foreach", "future", 
              "GGally", "ggplot2", "grid", "gridExtra", "Mass", 
              "panelView", "Rcpp")
install_all(packages)
library(fect)
data(fect)
ls()

##########################
# From: ./06-sens.Rmd
##########################


# install packages from CRAN
packages <- c("dplyr", "panelView", "ggplot2") # Removed HonestDiD, doParallel
install.packages(setdiff(packages, rownames(installed.packages())))  

# install most up-to-date "fect" from Github
if ("fect" %in% rownames(installed.packages()) == FALSE) {
  devtools:: install_github("xuyiqing/fect")
}

# install forked "HonestDiD" package compatible with "fect"
if ("HonestDiDFEct" %in% rownames(installed.packages()) == FALSE) {
  devtools:: install_github("lzy318/HonestDiDFEct") # This is used by fect_sens
}

library(dplyr)
library(fect)
library(panelView)
library(ggplot2)
library(HonestDiDFEct) # Required for fect_sens to work

data(fect)
data <- hh2019
head(data)

out.fect.placebo <- fect(nat_rate_ord~indirect, data = hh2019, 
                         index = c("bfs","year"),
                         method = 'fe', se = TRUE, 
                         placeboTest = TRUE, placebo.period = c(-2,0))

# Define post-treatment periods and sensitivity parameters for fect_sens
T.post <- 10 # Number of post-treatment periods based on original analysis
post_periods_vec <- 1:T.post

# Parameters for Relative Magnitude (RM) restriction
Mbar_vec_avg_rm <- seq(0, 1, by = 0.1)    # For average ATT plot
Mbar_vec_period_rm <- c(0, 0.5)          # For period-by-period ATT plot

# Parameters for Smoothness restriction
M_vec_avg_smooth <- seq(0, 0.25, by = 0.05) # For average ATT plot
M_vec_period_smooth <- c(0, 0.1)           # For period-by-period ATT plot

# Run sensitivity analysis using fect_sens
# This function augments out.fect.placebo with sensitivity results
out.fect.placebo <- fect_sens(
  fect.out      = out.fect.placebo,
  post.periods  = post_periods_vec,
  Mbarvec       = Mbar_vec_avg_rm,
  periodMbarvec = Mbar_vec_period_rm,
  Mvec          = M_vec_avg_smooth,
  periodMvec    = M_vec_period_smooth,
  parallel      = TRUE # Set to TRUE for parallel processing if desired
)

plot(out.fect.placebo,
     type = "sens",
     restrict = "rm",
     main = "Relative Magnitude Restriction")

plot(out.fect.placebo,
    type = "sens_es",
    restrict = "rm",
    main = "ATTs with Robust Confidence Sets (RM)",
    ylab = "Coefficients and 95% CI",
    xlim = c(-12,10), 
    ylim = c(-6,8), 
    show.count = TRUE)

plot(out.fect.placebo,
    type = "sens_es",
    restrict = "rm",
    main = "ATTs with Robust Confidence Sets (RM)",
    ylab = "Coefficients and 95% CI",
    xlim = c(-12,10), 
    ylim = c(-6,8), 
    show.count = TRUE,
    sens.colors = c("blue", "red"))

plot(out.fect.placebo,
    type = "sens",
    restrict = "sm",
    main = "Smoothness Restriction")

plot(out.fect.placebo,
    type = "sens_es",
    restrict = "sm",
    main = "ATTs with Robust Confidence Sets (Smoothness)",
    ylab = "Coefficients and 95% CI",
    xlim = c(-12,10), # Adjusted to match original detailed plot
    ylim = c(-12,15),
    show.count = TRUE)
