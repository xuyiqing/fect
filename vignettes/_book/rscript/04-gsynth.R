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
# From: ./04-gsynth.Rmd
##########################

set.seed(1234)
rm(list = ls())

library(fect)
data(fect)
ls()

head(simgsynth)

library(panelView)
panelview(Y ~ D, data = simgsynth,  index = c("id","time"), pre.post = TRUE) 

panelview(Y ~ D, data = simgsynth,  index = c("id","time"), type = "outcome") 

system.time(
out <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"), 
            method = "gsynth", force = "two-way", CV = TRUE, r = c(0, 5), 
            se = TRUE, nboots = 1000, vartype = 'parametric', 
            parallel = FALSE))

print(out)
out$est.att
out$est.avg
out$est.beta

system.time(
out <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"), method = "gsynth", force = "two-way", CV = TRUE, r = c(0, 5), se = TRUE, nboots = 1000,vartype = 'parametric', parallel = TRUE, cores = 16)
)

out2 <- fect(Y ~ D + X1 + X2, data = simgsynth,  index = c("id","time"), 
               method = "gsynth", force = "two-way", 
               CV = TRUE, r = c(0, 5), se = TRUE,
               vartype = "jackknife", 
               parallel = TRUE, cores = 16)


a <- plot(out) # by default, type = "gap"
print(a)

plot(out, theme.bw = FALSE) 

plot(out, connected = TRUE) 

plot(out, connected = TRUE, show.points = FALSE) 

plot(out, type = "gap", ylim = c(-3,12), xlab = "Period", 
     main = "Estimated ATT (FEct)")

plot(out, type = "status", yticklabels="0", 
     xticklabels=c("5", "10", "15","20", "25", "30") )

plot(out, type = "box", xlab = "time",
     xticklabels=c("-19", "-15", "-10", "-5","0","5","10") )

plot(out, type = "box", xlim = c(-15, 10), 
     xticklabels=c( "-15", "-10", "-5","0","5","10"))

plot(out,type = "calendar")

plot(out, type = "loadings")

plot(out, type = "factors", xlab = "Time")

plot(out, type = "counterfactual")

plot(out, type = "counterfactual", raw = "all")

plot(out, type = "counterfactual", raw = "band")

plot(out, type = "equiv", ylim = c(-5, 5))

plot(out, type = "equiv", show.stats =  FALSE)

plot(out, type = "equiv", stats.pos = c(-19, 4.5), ylim = c(-5, 5))

out.ife <- fect(Y ~ D + X1 + X2, data = simgsynth, index = c("id","time"), 
          force = "two-way", method = "ife", CV = TRUE, r = c(0, 5), 
          se = TRUE, nboots = 200, parallel = TRUE)

plot(out.ife, main = "Estimated ATT (EM)")

out.mc <- fect(Y ~ D + X1 + X2, data = simgsynth, 
               index = c("id","time"), 
                force = "two-way", method = "mc", CV = TRUE, 
                se = TRUE, nboots = 200, parallel = TRUE)


plot(out.mc, main = "Estimated ATT (MC)")

panelview(turnout ~ policy_edr, data = turnout, 
          index = c("abb","year"), pre.post = TRUE, 
          by.timing = TRUE) 

panelview(turnout ~ policy_edr, data = turnout, 
          index = c("abb","year"), type = "outcome", 
          main = "EDR Reform and Turnout", 
          by.group = TRUE)

out0 <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor, 
               data = turnout, index = c("abb","year"), 
               se = TRUE, method = "gsynth",
               r = 0, CV = FALSE, force = "two-way", 
               nboots = 1000, seed = 02139)

plot(out0, type = "gap", xlim = c(-15, 15))

out_turnout <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor, 
              data = turnout,  index = c("abb","year"), 
              se = TRUE, method = "gsynth", vartype = "parametric",
              r = c(0, 5), CV = TRUE, force = "two-way", 
              nboots = 1000, seed = 02139, keep.sims = TRUE)

dim(out_turnout$wgt.implied)
sort(out_turnout$wgt.implied[,8])

plot(out_turnout, xlim = c(-10, 5), ylim=c(-15,15))

plot(out_turnout, type = "status",xlab = "Year", ylab = "State", main = "Treatment Status", 
     xticklabels=c(1920, 1928, 1936, 1944, 1952, 1960, 
                   1968, 1976, 1984, 1992, 2000, 2008), xangle=10)

plot(out_turnout, type = "counterfactual")

plot(out_turnout, type = "counterfactual", id = "WI", main = "Wisconsin")

plot(out_turnout, type = "box", 
     xticklabels=c("-20", "-15", "-10", "-5","0","5","10"))

plot(out_turnout, type = "calendar", ylim = c(-15,15))

plot(out_turnout, type = "factors", xlab = "Year")

plot(out_turnout, type = "loadings")

set.seed(123456)
turnout.ub <- turnout[-c(which(turnout$abb=="WY")[1:15], 
                         sample(1:nrow(turnout),50,replace=FALSE)),]

panelview(turnout ~ policy_edr + policy_mail_in + policy_motor, 
          data = turnout.ub,  index = c("abb","year"), 
          pre.post = TRUE) 

out_ub <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor, 
              data = turnout.ub,  index = c("abb","year"), 
              se = TRUE, method = "gsynth", 
              r = c(0, 5), CV = TRUE, force = "two-way", 
              parallel = TRUE, min.T0 = 8, 
              nboots = 1000, seed = 02139)

plot(out_ub, type = "status",
     xticklabels=c(1920, 1928, 1936, 1944, 1952, 1960, 
                   1968, 1976, 1984, 1992, 2000, 2008),
     xangle=10)

plot(out_ub, type = "status", xlab = "Year", ylab = "State",
     main = "Treatment Status", id = out_ub$id[out_ub$tr],
     xlim = c(1920,2012), 
     xticklabels=c(1920, 1928, 1936, 1944, 1952, 1960,
                   1968, 1976, 1984, 1992, 2000, 2008))


plot(out_ub, type = "gap", ylim = c(-10, 20))

out.mc_ub <- fect(turnout ~ policy_edr + policy_mail_in + policy_motor, 
                 min.T0 = 8, data = turnout.ub,  
                 index = c("abb","year"), method = "mc", 
                 se = TRUE, nboots = 1000, seed = 02139)

plot(out.mc_ub, main = "Estimated ATT (MC)", ylim = c(-10, 20))