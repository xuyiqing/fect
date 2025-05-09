# Load required packages
library(HonestDiDFEct)
library(ggplot2)
# (Optional) If you have fect_sens_anlys defined in a separate script,
# you can source it here. Otherwise, make sure it's in your environment.
# source("path/to/fect_sens_anlys.R")

# Load HH2019 data
data(fect)


out.fect.placebo <- fect(nat_rate_ord~indirect, data = hh2019,
                         index = c("bfs","year"),
                         method = 'fe', se = TRUE,
                         placeboTest = TRUE, placebo.period = c(-2,0))

T.post <- 10 # 3 placebo periods, 10 post treatment periods
index.use <- which(rownames(out.fect.placebo$est.att) %in%
                     as.character(c(-2:T.post)))
beta.hat.p<- out.fect.placebo$est.att[index.use,1]
vcov.hat.p <- out.fect.placebo$att.vcov[index.use,index.use]
count <- out.fect.placebo$count[which(rownames(out.fect.placebo$est.att) %in%
                                        as.character(c(1:T.post)))]


honest.result.p <- HonestDiDFEct::createSensitivityResults_relativeMagnitudes(betahat = beta.hat.p,
                                                                              sigma = vcov.hat.p,
                                                                              numPrePeriods = 3,

                                                                              numPostPeriods = T.post,
                                                                              l_vec = count/sum(count),
                                                                              Mbarvec = seq(0,1,by=0.1))
#> Warning in .ARP_computeCI(betahat = betahat, sigma = sigma, numPrePeriods =
#> numPrePeriods, : CI is open at one of the endpoints; CI length may not be
#> accurate
#> Warning in .ARP_computeCI(betahat = betahat, sigma = sigma, numPrePeriods =
#> numPrePeriods, : CI is open at one of the endpoints; CI length may not be
#> accurate

originalResults.p <- HonestDiDFEct::constructOriginalCS(betahat = beta.hat.p,
                                                        sigma = vcov.hat.p,

                                                        numPrePeriods = 3,
                                                        numPostPeriods = T.post,
                                                        l_vec = count/sum(count))

HonestDiDFEct::createSensitivityPlot_relativeMagnitudes(honest.result.p, originalResults.p)


dte_base <- rep(0,T.post)
dte_output <- cbind.data.frame(lb = c(),
                               ub = c(),
                               method = c(),
                               Delta = c(),
                               Mbar = c(),
                               post = c())

for(t.post in c(1:T.post)){
  dte_l <- dte_base
  dte_l[t.post] <- 1
  honest.dte <- HonestDiDFEct::createSensitivityResults_relativeMagnitudes(betahat = beta.hat.p,
                                                                           sigma = vcov.hat.p,
                                                                           numPrePeriods = 3,

                                                                           numPostPeriods = T.post,
                                                                           l_vec = dte_l,
                                                                           Mbarvec = c(0,0.5))

  honest.dte <- as.data.frame(honest.dte)
  honest.dte$post <- t.post
  dte_output <- rbind(dte_output,honest.dte)
}
# event study plot
fect.output.p <- as.data.frame(out.fect.placebo$est.att)
fect.output.p$Time <- as.numeric(rownames(fect.output.p))
p.placebo.honest <- esplot(fect.output.p,Period = 'Time',Estimate = 'ATT',
                           SE = 'S.E.',CI.lower = "CI.lower",
                           main = "ATTs with Robust Confidence Sets (RM)",
                           ylab = "Coefficients and 95% CI",
                           CI.upper = "CI.upper",
                           xlim = c(-12,10),
                           ylim = c(-6,8),
                           Count = 'count', show.count = T)
dte_output <- as.data.frame(dte_output)

# add confidence sets
p.placebo.honest <- p.placebo.honest +
  geom_linerange(aes(x=post+0.2,ymin=lb,ymax=ub),
                 data = dte_output[which(dte_output$Mbar==0.5),],
                 color = 'maroon1',linewidth = 1)

# estimates assuming no PT violation
p.placebo.honest <- p.placebo.honest +
  geom_linerange(aes(x=post+0.2,ymin=lb,ymax=ub),
                 data = dte_output[which(dte_output$Mbar==0),],
                 color = '#228B22',linewidth = 1)

# placebo estimates
p.placebo.honest <- p.placebo.honest +
  geom_linerange(aes(x=Time,ymin=CI.lower,ymax=CI.upper),
                 data = fect.output.p[which(fect.output.p$Time%in%c(-2:0)),],
                 color = 'blue',linewidth = 1)
p.placebo.honest




dte_base <- rep(0,T.post)
dte_output <- cbind.data.frame(lb = c(),
                               ub = c(),
                               method = c(),
                               Delta = c(),
                               Mbar = c(),
                               post = c())

for(t.post in c(1:T.post)){
  dte_l <- dte_base
  dte_l[t.post] <- 1
  honest.dte <- HonestDiDFEct::createSensitivityResults(betahat = beta.hat.p,
                                                        sigma = vcov.hat.p,method = 'C-LF',
                                                        numPrePeriods = 3,
                                                        numPostPeriods = T.post,
                                                        l_vec = dte_l,
                                                        Mvec = c(0,0.1))

  honest.dte <- as.data.frame(honest.dte)
  honest.dte$post <- t.post
  dte_output <- rbind(dte_output,honest.dte)
}


# prepare for plotting
fect.output.p <- as.data.frame(out.fect.placebo$est.att)
fect.output.p$Time <- as.numeric(rownames(fect.output.p))
p.placebo.honest <- esplot(fect.output.p,Period = 'Time',Estimate = 'ATT',
                           SE = 'S.E.',CI.lower = "CI.lower",
                           main = "ATTs with Robust Confidence Sets (Smoothness)",
                           ylab = "Coefficients and 95% CI",
                           CI.upper = "CI.upper",xlim = c(-12,10),
                           ylim = c(-15,15),Count = 'count',show.count = T)
dte_output <- as.data.frame(dte_output)

# estimates assuming potential PT violation w/ a linear trend + partial deviation
p.placebo.honest <- p.placebo.honest +
  geom_linerange(aes(x=post+0.2,ymin=lb,ymax=ub),
                 data = dte_output[which(dte_output$M==0.1),],
                 color = 'skyblue',linewidth = 1)

# estimates assuming a strictly linear PT violation
p.placebo.honest <- p.placebo.honest +
  geom_linerange(aes(x=post+0.2,ymin=lb,ymax=ub),
                 data = dte_output[which(dte_output$M==0),],
                 color = 'orange',linewidth = 1)

# estimates assuming no PT violation
p.placebo.honest <- p.placebo.honest +
  geom_linerange(aes(x=Time,ymin=CI.lower,ymax=CI.upper),
                 data = fect.output.p[which(fect.output.p$Time>0),],
                 color = 'black',linewidth = 1)

# placebo estimates
p.placebo.honest <- p.placebo.honest +
  geom_linerange(aes(x=Time,ymin=CI.lower,ymax=CI.upper),
                 data = fect.output.p[which(fect.output.p$Time%in%c(-2:0)),],
                 color = 'blue',linewidth = 1) + xlim(-12.5,10.5)
#> Scale for x is already present.
#> Adding another scale for x, which will replace the existing scale.
p.placebo.honest
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_segment()`).
