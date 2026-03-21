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
# From: ./05-panel.Rmd
##########################

# install packages from CRAN
packages <- c("dplyr", "fixest", "did", "didimputation", 
              "panelView", "ggplot2", "bacondecomp", "HonestDiD",
              "DIDmultiplegtDYN", "PanelMatch")
install.packages(setdiff(packages, rownames(installed.packages())))  

# install most up-to-date "fect" from Github
if ("fect" %in% rownames(installed.packages()) == FALSE) {
  devtools:: install_github("xuyiqing/fect")
}


library(dplyr)
library(readstata13)
library(fixest)
library(did)
library(fect)
library(panelView)
library(PanelMatch)
library(ggplot2)
library(bacondecomp)
library(fect)
library(didimputation)
library(doParallel)
library(HonestDiD)
library(HonestDiDFEct)
library(DIDmultiplegtDYN) # may require XQuartz 

data(fect)
data <- hh2019
head(data)

panelview(nat_rate_ord ~ indirect, data = data, index = c("bfs","year"), 
  xlab = "Year", ylab = "Unit", display.all = T,
  gridOff = TRUE, by.timing = TRUE)

panelview(data = data,Y='nat_rate_ord',
          D='indirect',index=c("bfs","year"),
          by.timing = TRUE, display.all = TRUE,
          type = "outcome", by.cohort = TRUE)

# remember to cluster standard errors
model.twfe.0 <- feols(nat_rate_ord~indirect|bfs+year,
                      data=data, cluster = "bfs") 
print(model.twfe.0)

data.complete <- data[which(!is.na(data$nat_rate_ord)),] # bacon requires no missingness in the data
df_bacon <- bacon(nat_rate_ord~indirect,
                  data = data.complete,
                  id_var = "bfs",
                  time_var = "year")
ggplot(df_bacon) +
   aes(x = weight, y = estimate, shape = factor(type), color = factor(type)) +
   labs(x = "Weight", y = "Estimate", shape = "Type", color = 'Type') +
   geom_point()

print(aggregate(df_bacon$estimate * df_bacon$weight, 
                list(df_bacon$type), FUN=sum))

# drop always treated units
df <- as.data.frame(data %>% 
                      group_by(bfs) %>%
                      mutate(treatment_mean = mean(indirect,na.rm = TRUE)))
df.use <- df[which(df$treatment_mean<1),]

# Re-estimate TWFE on this Sub-sample
model.twfe.1 <- feols(nat_rate_ord~indirect|bfs+year,
                      data=df.use, cluster = "bfs")
print(model.twfe.1)

df.use <- get.cohort(df.use, D = "indirect", index=c("bfs","year"), 
                     start0 = TRUE)
head(df.use[,-5],19)

# Dynamic TWFE
df.twfe <- df.use
# drop always treated units
df.twfe$treat <- as.numeric(df.twfe$treatment_mean>0) 
df.twfe[which(is.na(df.twfe$Time_to_Treatment)),'Time_to_Treatment'] <- 0 # can be an arbitrary value
twfe.est <- feols(nat_rate_ord ~ i(Time_to_Treatment, treat, ref = -1)| bfs + year, 
                  data = df.twfe, cluster = "bfs")
twfe.output <- as.matrix(twfe.est$coeftable)
print(round(twfe.output, 3))

twfe.output <- as.data.frame(twfe.output)
twfe.output$Time <- c(c(-18:-2),c(0:17))+1 
p.twfe <- esplot(twfe.output,Period = 'Time',Estimate = 'Estimate',
                               SE = 'Std. Error', xlim = c(-12,10))
p.twfe

twfe.output <- as.data.frame(twfe.est$coeftable)
twfe.output$Time <- c(c(-18:-2),c(0:17)) 
p.twfe <- esplot(twfe.output, Period = 'Time',
                 Estimate = 'Estimate', SE = 'Std. Error', 
                 xlim = c(-12,10),start0 = TRUE)
p.twfe

df.st <- NULL
target.cohorts <- setdiff(unique(df.use$Cohort),"Control")
k <- 1
for(cohort in target.cohorts){
  df.sub <- df.use[which(df.use$Cohort%in%c(cohort,"Control")),]
  df.sub$stack <- k
  df.st <- rbind(df.st,df.sub)
  k <- k + 1
}
df.st$st_unit <- as.numeric(factor(paste0(df.st$stack,'-',df.st$bfs)))
df.st$st_year <- as.numeric(factor(paste0(df.st$stack,'-',df.st$year)))
model.st <- feols(nat_rate_ord~indirect|st_unit+st_year,
                  data=df.st, cluster = "st_unit")

print(model.st)

df.st$treat <- as.numeric(df.st$treatment_mean>0)
df.st[which(is.na(df.st$Time_to_Treatment)),'Time_to_Treatment'] <- 1000 
# note, this "1000" can be arbitrary value

st.est <- feols(nat_rate_ord ~ 
                  i(Time_to_Treatment, treat, ref = -1)| st_unit + 
                  st_year,data = df.st,cluster = "st_unit")

# make plot
st.output <- as.data.frame(st.est$coeftable)
st.output$Time <- c(c(-18:-2),c(0:17))+1 
p.st <- esplot(st.output,Period = 'Time',Estimate = 'Estimate',
                               SE = 'Std. Error', xlim = c(-12,10))
p.st

knitr::include_graphics("fig/fig_iw.png")

df.sa <- df.use
df.sa[which(is.na(df.sa$FirstTreat)),"FirstTreat"] <- 1000 
# above, replace NA with an arbitrary number 

model.sa.1 <- feols(nat_rate_ord~sunab(FirstTreat,year)|bfs+year,
                    data = df.sa, cluster = "bfs")
summary(model.sa.1,agg = "ATT")

sa.output <- as.data.frame(as.matrix(model.sa.1$coeftable))
sa.output$Time <- c(c(-18:-2),c(0:17)) + 1
p.sa <- esplot(sa.output,Period = 'Time',Estimate = 'Estimate',
                             SE = 'Std. Error', xlim = c(-12,10))
p.sa

df.cs <- df.use
df.cs[which(is.na(df.cs$FirstTreat)),"FirstTreat"] <- 0 # replace NA with 0
cs.est.1 <- att_gt(yname = "nat_rate_ord",
                 gname = "FirstTreat",
                 idname = "bfs",
                 tname = "year",
                 xformla = ~1,
                 control_group = "nevertreated",
                 allow_unbalanced_panel = TRUE,
                 data = df.cs,
                 est_method = "reg")
cs.est.att.1 <- aggte(cs.est.1, type = "simple", na.rm=T, bstrap = F)
print(cs.est.att.1)

cs.att.1 <- aggte(cs.est.1, type = "dynamic",
                  bstrap=FALSE, cband=FALSE, na.rm=T) 
print(cs.att.1)

cs.output <- cbind.data.frame(Estimate = cs.att.1$att.egt,
                              SE = cs.att.1$se.egt,
                              time = cs.att.1$egt + 1)
p.cs.1 <- esplot(cs.output,Period = 'time',Estimate = 'Estimate',
                               SE = 'SE', xlim = c(-12,10))
p.cs.1

cs.est.1.u <- att_gt(yname = "nat_rate_ord",
                 gname = "FirstTreat",
                 idname = "bfs",
                 tname = "year",
                 xformla = ~1,
                 control_group = "nevertreated",
                 allow_unbalanced_panel = TRUE,
                 data = df.cs,
                 est_method = "reg", 
                 base_period = "universal")
cs.att.1.u <- aggte(cs.est.1.u, type = "dynamic",
                    bstrap=FALSE, cband=FALSE, na.rm=T) 
cs.output.u <- cbind.data.frame(Estimate = cs.att.1.u$att.egt,
                                SE = cs.att.1.u$se.egt,
                                time = cs.att.1.u$egt + 1)
p.cs.1.u <- esplot(cs.output.u,Period = 'time',Estimate = 'Estimate',
                               SE = 'SE', xlim = c(-12,10))
p.cs.1.u

knitr::include_graphics("fig/fig_cs.png")

cs.est.2 <- att_gt(yname = "nat_rate_ord",
                   gname = "FirstTreat",
                   idname = "bfs",
                   tname = "year",
                   xformla = ~1,
                   control_group = "notyettreated",
                   allow_unbalanced_panel = TRUE,
                   data = df.cs,
                   est_method = "reg")
cs.est.att.2 <- aggte(cs.est.2, type = "simple",na.rm=T, bstrap = F)
print(cs.est.att.2)

cs.est.2.u <- att_gt(yname = "nat_rate_ord", gname = "FirstTreat",
                     idname = "bfs", tname = "year", xformla = ~1,
                     control_group = "notyettreated",
                     allow_unbalanced_panel = TRUE,
                     data = df.cs, est_method = "reg", 
                     base_period = "universal")

cs.att.2.u <- aggte(cs.est.2.u, type = "dynamic",
                    bstrap=FALSE, cband=FALSE, na.rm=T) 

# plot
cs.output.u <- cbind.data.frame(Estimate = cs.att.2.u$att.egt,
                                SE = cs.att.2.u$se.egt,
                                time = cs.att.2.u$egt + 1)
p.cs.2.u <- esplot(cs.output.u,Period = 'time',Estimate = 'Estimate',
                               SE = 'SE', xlim = c(-12,10))
p.cs.2.u

knitr::include_graphics("fig/fig_pm.png")

didm.results <- did_multiplegt_dyn(
      df = df.use,
      outcome = "nat_rate_ord",
      group = "bfs",
      controls = NULL,
      time = "year",
      treatment = "indirect",
      effects = 12, 
      placebo = 9,
      cluster = "bfs",
      graph_off = TRUE
    )
print(didm.results)

T.post <- dim(didm.results$results$Effects)[1]
T.pre <- dim(didm.results$results$Placebos)[1]
didm.vis <- rbind(didm.results$results$Placebos,didm.results$results$Effects)
didm.vis <- as.data.frame(didm.vis)
didm.vis[,'Time'] <- c(c(-1:-(T.pre)),c(1:T.post))
est.dynamic <- didm.vis[,c(9,1,2,3,4)]
colnames(est.dynamic) <- c("T","estimate","se","lb","ub")
p.didm <- esplot(est.dynamic,Period = 'T',Estimate = 'estimate',
                               SE = 'se', xlim = c(-9, 9))
p.didm

df.pm <- df.use
# we need to convert the unit and time indicator to integer
df.pm[,"bfs"] <- as.integer(as.factor(df.pm[,"bfs"]))
df.pm[,"year"] <- as.integer(as.factor(df.pm[,"year"]))
df.pm <- df.pm[,c("bfs","year","nat_rate_ord","indirect")]

# Pre-processes and balances panel data
df.pm <- PanelData(panel.data = df.pm,
                    unit.id = "bfs",
                    time.id = "year",
                    treatment = "indirect",
                    outcome = "nat_rate_ord")

PM.results <- PanelMatch(lag=3, 
                         refinement.method = "none", 
                         panel.data = df.pm, 
                         qoi = "att", 
                         lead = c(0:3), 
                         match.missing = TRUE)

## For pre-treatment dynamic effects
PM.results.placebo <- PanelMatch(lag=3, 
                         refinement.method = "none", 
                         panel.data = df.pm, 
                         qoi = "att", 
                         lead = c(0:3), 
                         match.missing = TRUE,
                         placebo.test = TRUE)


# ATT
PE.results.pool <- PanelEstimate(PM.results, panel.data = df.pm, pooled = TRUE)
summary(PE.results.pool)

# Dynamic Treatment Effects
PE.results <- PanelEstimate(PM.results, panel.data = df.pm)
PE.results.placebo <- placebo_test(PM.results.placebo, panel.data = df.pm, plot = F)

# obtain lead and lag (placebo) estimates
est_lead <- as.vector(PE.results$estimate)
est_lag <- as.vector(PE.results.placebo$estimates)
sd_lead <- apply(PE.results$bootstrapped.estimates,2,sd)
sd_lag <- apply(PE.results.placebo$bootstrapped.estimates,2,sd)
coef <- c(est_lag, 0, est_lead)
sd <- c(sd_lag, 0, sd_lead)
pm.output <- cbind.data.frame(ATT=coef, se=sd, t=c(-2:4))

# plot
p.pm <- esplot(data = pm.output,Period = 't',
               Estimate = 'ATT',SE = 'se')
p.pm

knitr::include_graphics("fig/fig_fect.png")

out.fect <- fect(nat_rate_ord~indirect, data = df, 
                 index = c("bfs","year"),
                 method = 'fe', se = TRUE)
print(out.fect$est.avg)

fect.output <- as.matrix(out.fect$est.att)
head(fect.output)

df.impute <- df.use
df.impute[which(is.na(df.impute$FirstTreat)),"FirstTreat"] <- 0 
# above, replace NA with 0

out.impute <- did_imputation(data = df.impute,
                               yname = "nat_rate_ord",
                               gname = "FirstTreat",
                               tname = "year",
                               idname = "bfs",
                               cluster_var = "bfs")
out.impute

fect.output <- as.data.frame(fect.output)
fect.output$Time <- c(-17:18)
p.fect <- esplot(fect.output,Period = 'Time',Estimate = 'ATT',
                   SE = 'S.E.',CI.lower = "CI.lower", 
                   CI.upper = 'CI.upper',xlim = c(-12,10))
p.fect

model.impute <- did_imputation(data = df.impute,
                               yname = "nat_rate_ord",
                               gname = "FirstTreat",
                               tname = "year",
                               idname = "bfs",
                               cluster_var = "bfs",
                               pretrends = c(-13:-1),
                               horizon = TRUE)
model.impute$term <- as.numeric(model.impute$term)+1 
# above, set 1 as the first post-treatment period

# plot
to_plot <- as.data.frame(model.impute)
esplot(data=to_plot,Period = "term", 
       Estimate = 'estimate', SE = 'std.error',
       xlim = c(-12,10))
out.impute

out.fect.balance <- fect(nat_rate_ord~indirect, data = df, 
                         index = c("bfs","year"),
                         method = 'fe', se = TRUE, 
                         balance.period = c(-2,4))
# att
print(out.fect.balance$est.balance.avg)

# event study plot
fect.balance.output <- as.data.frame(out.fect.balance$est.balance.att)
fect.balance.output$Time <- c(-2:4)
p.fect.balance <- esplot(fect.balance.output,Period = 'Time',
                         Estimate = 'ATT', SE = 'S.E.',
                         CI.lower = "CI.lower", 
                         CI.upper = 'CI.upper')
p.fect.balance

data(fect)
data <- gs2020
data$cycle <- as.integer(as.numeric(data$cycle/2))
head(data)

y <- "general_sharetotal_A_all"
d <- "cand_A_all"
unit <- "district_final"
time <- "cycle"
controls <- c("cand_H_all", "cand_B_all")
index <- c("district_final", "cycle")

panelview(Y=y, D=d, X=controls, index = index, data = data, 
          xlab = "Time Period", ylab = "Unit", gridOff = TRUE, 
          by.timing = TRUE, cex.legend=5, cex.axis= 5, 
          cex.main = 10, cex.lab = 5)


model.twfe <- feols(general_sharetotal_A_all ~ cand_A_all + 
                      cand_H_all + cand_B_all | district_final + cycle,
                    data=data, cluster = "district_final") 
summary(model.twfe)

data_cohort <- get.cohort(data, index = index, D=d,start0 = TRUE)
# Generate a dummy variable treat
data_cohort$treat <- 0
data_cohort[which(data_cohort$Cohort!='Control'),'treat'] <- 1
data_cohort[which(is.na(data_cohort$Time_to_Treatment)), "treat"] <- 0

# remove observations that starts with treated status
remove <- intersect(which(is.na(data_cohort$Time_to_Treatment)),
                    which(data_cohort[,d]==1)) 
if(length(remove)>0){data_cohort <- data_cohort[-remove,]}

# replace missingness in Time_to_Treatment with an arbitrary number
data_cohort[which(is.na(data_cohort$Time_to_Treatment)), "Time_to_Treatment"] <- 999 

twfe.est <- feols(general_sharetotal_A_all ~ 
                    i(Time_to_Treatment, treat, ref = -1) + 
                    cand_H_all +cand_B_all | district_final + cycle,  
                  data = data_cohort, cluster = "district_final")

twfe.output <- as.data.frame(twfe.est$coeftable[c(1:25),])
twfe.output$Time <- c(c(-16:-2),c(0:9)) + 1 

# plot
p.twfe <- esplot(twfe.output,Period = 'Time',Estimate = 'Estimate',
                               SE = 'Std. Error', xlim = c(-15,1))
p.twfe

df.pm <- data_cohort
# we need to convert the unit and time indicator to integer
df.pm[,"district_final"] <- as.integer(as.factor(df.pm[,"district_final"]))
df.pm[,"cycle"] <- as.integer(as.factor(df.pm[,"cycle"]))
df.pm <- df.pm[,c("district_final","cycle","cand_A_all", 
                  "general_sharetotal_A_all")]

# Pre-processes and balances panel data
df.pm <- PanelData(panel.data = df.pm,
                    unit.id = "district_final",
                    time.id = "cycle",
                    treatment = "cand_A_all",
                    outcome = "general_sharetotal_A_all")

PM.results <- PanelMatch(lag=4, 
                         refinement.method = "none", 
                         panel.data = df.pm, 
                         qoi = "att", 
                         lead = 0, 
                         match.missing = TRUE)

## For pre-treatment dynamic effects
PM.results.placebo <- PanelMatch(lag=4, 
                         refinement.method = "none", 
                         panel.data = df.pm, 
                         qoi = "att", 
                         lead = 0, 
                         match.missing = TRUE,
                         placebo.test = TRUE)


PE.results.pool <- PanelEstimate(PM.results, panel.data = df.pm, pooled = TRUE)
summary(PE.results.pool)

# Dynamic Treatment Effects
PE.results <- PanelEstimate(PM.results, panel.data = df.pm)
PE.results.placebo <- placebo_test(PM.results.placebo, panel.data = df.pm,
                                   plot = FALSE)

est_lead <- as.vector(PE.results$estimate)
est_lag <- as.vector(PE.results.placebo$estimates)
sd_lead <- apply(PE.results$bootstrapped.estimates,2,sd)
sd_lag <- apply(PE.results.placebo$bootstrapped.estimates,2,sd)
coef <- c(est_lag, 0, est_lead)
sd <- c(sd_lag, 0, sd_lead)
pm.output <- cbind.data.frame(ATT=coef, se=sd, t=c(-3:1))

# plot
p.pm <- esplot(data = pm.output,Period = 't',
               Estimate = 'ATT',SE = 'se')
p.pm

model.fect <- fect(Y = "general_sharetotal_A_all", D = "cand_A_all", 
                   X= c("cand_H_all", "cand_B_all"), data = data, 
                   method = "fe", index = index, se = TRUE, 
                   parallel = TRUE, seed = 1234, force = "two-way")

print(model.fect$est.avg)

fect.output <- as.data.frame(model.fect$est.att)
fect.output$Time <- c(-15:10)
p.fect <- esplot(fect.output,Period = 'Time',Estimate = 'ATT',
                   SE = 'S.E.',CI.lower = "CI.lower", 
                   CI.upper = 'CI.upper', xlim = c(-15,1))
p.fect

plot(model.fect)

plot(model.fect, type = 'exit')

out.fect.p <- fect(Y = y, X = controls, D = d, data = data, index = index,
                   method = 'fe', se = TRUE, placeboTest = TRUE,
                   placebo.period = c(-2,0))

plot(out.fect.p, proportion = 0.1, stats = "placebo.p")

out.fect.c <- fect(Y = y, X = controls, D = d, data = data, index = index,
                   method = 'fe', se = TRUE, carryoverTest = TRUE, carryover.period = c(1,2))

# plot
plot(out.fect.c,  stats = "carryover.p", ylim = c(-0.15, 0.20))

out.fect.balance <- fect(Y = y, X = controls, D = d, data = data, 
                         index = index, method = 'fe', se = TRUE,
                         balance.period = c(-3,1))

# att
print(out.fect.balance$est.balance.avg)

# event study plot
fect.balance.output <- as.data.frame(out.fect.balance$est.balance.att)
fect.balance.output$Time <- c(-3:1)
p.fect.balance <- esplot(fect.balance.output,Period = 'Time',Estimate = 'ATT',
                   SE = 'S.E.',CI.lower = "CI.lower", 
                   CI.upper = 'CI.upper')
p.fect.balance

res_st <- did_wrapper(
  data   = hh2019,
  Y      = "nat_rate_ord",
  D      = "indirect",
  index  = c("bfs", "year"),
  method = "st",
  se     = "default"
)
print(res_st)

res_st <- did_wrapper(
  data   = hh2019,
  Y      = "nat_rate_ord",
  D      = "indirect",
  index  = c("bfs", "year"),
  method = "st",
  se     = "boot"
)
print(res_st)

esplot(data = res_st, main = "Stacked DID", xlim = c(-12,10))
