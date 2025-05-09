# Load required packages
library(HonestDiDFEct)

# (Optional) If you have fect_sens_anlys defined in a separate script,
# you can source it here. Otherwise, make sure it's in your environment.
# source("path/to/fect_sens_anlys.R")

# Load HH2019 data
data(fect)

out.fect <- fect(
  nat_rate_ord ~ indirect,
  data  = hh2019,
  index = c("bfs", "year"),
  method        = 'fe',
  se            = TRUE,
)

# # Fit fect with placeboTest = TRUE
# out.fect.placebo3 <- fect(
#   nat_rate_ord ~ indirect,
#   data  = hh2019,
#   index = c("bfs", "year"),
#   method        = 'fe',
#   se            = TRUE,
#   placeboTest   = TRUE,
#   placebo.period = c(-2, 0)
# )

out.fect.placebo.ife <- fect(
  Y ~ D + X1 + X2,
  data  = simgsynth,
  index = c("id","time"),
  r= 2,
  force = "two-way",
  nboots = 1000,
  parallel =  TRUE,
  method        = 'ife',
  se            = TRUE,
  placeboTest   = TRUE,
  placebo.period = c(-2, 0)
  ,vartype = 'parametric'
)

out.fect.placebo <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"),
        method = "gsynth", force = "two-way", r= 2,
        se = TRUE, nboots = 1000,vartype = 'parametric', parallel = TRUE, cores = 16,  placeboTest   = TRUE,placebo.period = c(-2, 0))

# Before we run the sensitivity function, ensure that 'att.vcov' has dimnames
# to match 'est.att'. You may or may not need this, depending on your version.
# If needed:
if (is.null(dimnames(out.fect.placebo$att.vcov))) {
  dimnames(out.fect.placebo$att.vcov) <- list(
    rownames(out.fect.placebo$est.att),
    rownames(out.fect.placebo$est.att)
  )
}

# Run the sensitivity analysis (the updated function returns the fect object itself)
out.fect.placebo2 <- fect_sens(
  fect.out     = out.fect.placebo,
  post.periods = 1:10,
  periodMbarvec      = c(1),
  parallel = FALSE
)

out.fect.placebo.ife2 <- fect_sens(
  fect.out     = out.fect.placebo.ife,
  post.periods = 1:10,
  periodMbarvec      = c(1),
)



# out.fect.placebo2 <- fect_sens(
#   fect.out     = out.fect.placebo,
#   post.periods = 1:10,
#   Mbarvec      = seq(0,1, by=0.1),
#   periodMbarvec      = seq(0,1, by=0.1),
#   Mvec         = seq(0,0.25,0.05),
#   periodMvec         = seq(0,0.25,0.05)
# )

# -------------------------------------------------------------------
# Inspect the new sub-lists in out.fect.placebo
# -------------------------------------------------------------------


HonestDiDFEct::createSensitivityPlot_relativeMagnitudes( out.fect.placebo2$sensitivity.rm$results, out.fect.placebo2$sensitivity.rm$original)

# Weighted RM:
plot(out.fect.placebo2, type = "sens_rm")

# Weighted Smoothness:
plot(out.fect.placebo2, type = "sens_smooth")

# Period-by-Period RM:
plot(out.fect.placebo2, xlim = c(-12,10), type = "sens_smooth_gaps", theme = "grayscale")

# Period-by-Period Smoothness:
plot(out.fect.placebo2,  xlim = c(-12,10),type = "sens")

esplot(data = out.fect.placebo$est.att, connected = TRUE, highlight.periods = -2:0, highlight.colors = rep("blue",3), show.count = TRUE)

esplot(data = out.fect.placebo$est.att,  highlight.periods = -2:0, highlight.colors = rep("blue",3), show.count = TRUE)
plot(out.fect.placebo2)
