## ----.common, include = FALSE-------------------------------------------------
source("_common.R")


## ----install-packages, message = FALSE, warning = FALSE-----------------------
# install packages from CRAN
packages <- c("dplyr", "fixest", "did", "didimputation",
              "panelView", "ggplot2", "bacondecomp", "HonestDiD",
              "DIDmultiplegtDYN", "PanelMatch", "readstata13")
install.packages(setdiff(packages, rownames(installed.packages())))  

# install most up-to-date "fect" from Github
if ("fect" %in% rownames(installed.packages()) == FALSE) {
  devtools:: install_github("xuyiqing/fect")
}

# install forked "HonestDiD" package compatible with "fect"
if ("HonestDiDFEct" %in% rownames(installed.packages()) == FALSE) {
  devtools:: install_github("lzy318/HonestDiDFEct")
}


## ----load-libraries, message = FALSE, warning = FALSE-------------------------
library(dplyr)
library(readstata13)
library(fixest)
library(did)
library(panelView)
library(PanelMatch)
library(ggplot2)
library(bacondecomp)
library(didimputation)
library(doParallel)
library(HonestDiD)
library(HonestDiDFEct)
has_polars <- requireNamespace("polars", quietly = TRUE)
if (has_polars) {
  library(polars)
  library(DIDmultiplegtDYN) # requires polars; may require XQuartz for rgl
}


## ----load-hh2019, message = FALSE, warning = FALSE----------------------------
data(hh2019)
data <- hh2019
head(data)


## ----hh_panelview_treat, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
panelview(nat_rate_ord ~ indirect, data = data, index = c("bfs","year"), 
  xlab = "Year", ylab = "Unit", display.all = T,
  gridOff = TRUE, by.timing = TRUE)


## ----hh_panelview_cohort, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5,cache=TRUE----
panelview(data = data,Y='nat_rate_ord',
          D='indirect',index=c("bfs","year"),
          by.timing = TRUE, display.all = TRUE,
          type = "outcome", by.cohort = TRUE)


## ----hh_twfe1, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5----
# remember to cluster standard errors
model.twfe.0 <- feols(nat_rate_ord~indirect|bfs+year,
                      data=data, cluster = "bfs") 
print(model.twfe.0)


## ----hh_bacon, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
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


## ----hh_twfe2, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
# drop always treated units
df <- as.data.frame(data %>% 
                      group_by(bfs) %>%
                      mutate(treatment_mean = mean(indirect,na.rm = TRUE)))
df.use <- df[which(df$treatment_mean<1),]

# Re-estimate TWFE on this Sub-sample
model.twfe.1 <- feols(nat_rate_ord~indirect|bfs+year,
                      data=df.use, cluster = "bfs")
print(model.twfe.1)


## ----hh_cohort, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
df.use <- get.cohort(df.use, D = "indirect", index=c("bfs","year"), 
                     start0 = TRUE)
head(df.use[,-5],19)


## ----hh_twfeplot, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
# Dynamic TWFE
df.twfe <- df.use
# drop always treated units
df.twfe$treat <- as.numeric(df.twfe$treatment_mean>0) 
df.twfe[which(is.na(df.twfe$Time_to_Treatment)),'Time_to_Treatment'] <- 0 # can be an arbitrary value
twfe.est <- feols(nat_rate_ord ~ i(Time_to_Treatment, treat, ref = -1)| bfs + year, 
                  data = df.twfe, cluster = "bfs")
twfe.output <- as.matrix(twfe.est$coeftable)
print(round(twfe.output, 3))


## ----hh_twfeplot2, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
twfe.output <- as.data.frame(twfe.output)
twfe.output$Time <- c(c(-18:-2),c(0:17))+1 
p.twfe <- esplot(twfe.output,Period = 'Time',Estimate = 'Estimate',
                               SE = 'Std. Error', xlim = c(-12,10))
p.twfe


## ----hh_twfeplot3, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
twfe.output <- as.data.frame(twfe.est$coeftable)
twfe.output$Time <- c(c(-18:-2),c(0:17))
p.twfe <- esplot(twfe.output, Period = 'Time',
                 Estimate = 'Estimate', SE = 'Std. Error',
                 xlim = c(-12,10),start0 = TRUE)
p.twfe


## ----hh_twfeplot_connected, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
p.twfe.connected <- esplot(twfe.output, Period = 'Time',
                           Estimate = 'Estimate', SE = 'Std. Error',
                           xlim = c(-12, 10), start0 = TRUE,
                           connected = TRUE,
                           main = "TWFE event study (connected)")
p.twfe.connected


## ----hh_st, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
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


## ----hh_stplot, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
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


## ----fig-iw-illustration, echo=FALSE, out.width="50%", fig.align="center"-----
knitr::include_graphics("fig/fig_iw.png")


## ----hh_sa, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
df.sa <- df.use
df.sa[which(is.na(df.sa$FirstTreat)),"FirstTreat"] <- 1000 
# above, replace NA with an arbitrary number 

model.sa.1 <- feols(nat_rate_ord~sunab(FirstTreat,year)|bfs+year,
                    data = df.sa, cluster = "bfs")
summary(model.sa.1,agg = "ATT")


## ----hh_saplot, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
sa.output <- as.data.frame(as.matrix(model.sa.1$coeftable))
sa.output$Time <- c(c(-18:-2),c(0:17)) + 1
p.sa <- esplot(sa.output,Period = 'Time',Estimate = 'Estimate',
                             SE = 'Std. Error', xlim = c(-12,10))
p.sa


## ----hh_cs1, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
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


## ----hh_csplot1, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
cs.att.1 <- aggte(cs.est.1, type = "dynamic",
                  bstrap=FALSE, cband=FALSE, na.rm=T) 
print(cs.att.1)


## ----hh_csplot1a, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
cs.output <- cbind.data.frame(Estimate = cs.att.1$att.egt,
                              SE = cs.att.1$se.egt,
                              time = cs.att.1$egt + 1)
p.cs.1 <- esplot(cs.output,Period = 'time',Estimate = 'Estimate',
                               SE = 'SE', xlim = c(-12,10))
p.cs.1


## ----hh_csplot1b, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
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


## ----fig-cs-illustration, echo=FALSE, out.width="50%", fig.align="center"-----
knitr::include_graphics("fig/fig_cs.png")


## ----hh_cs2, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
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


## ----hh_csplot2b, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
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


## ----fig-pm-illustration, echo=FALSE, out.width="50%", fig.align="center"-----
knitr::include_graphics("fig/fig_pm.png")






## ----hh_pm, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
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



## ----hh_pm1, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
# ATT
PE.results.pool <- PanelEstimate(PM.results, panel.data = df.pm, pooled = TRUE)
summary(PE.results.pool)


## ----hh_pm2, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
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


## ----fig-fect-illustration, echo=FALSE, out.width="50%", fig.align="center"----
knitr::include_graphics("fig/fig_fect.png")


## ----hh_fect, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
out.fect <- fect(nat_rate_ord~indirect, data = df, 
                 index = c("bfs","year"),
                 method = 'fe', se = TRUE)
print(out.fect$est.avg)


## ----hh_fectplot, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
fect.output <- as.matrix(out.fect$est.att)
head(fect.output)


## ----hh_impute, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
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


## ----hh_fectplot2, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
fect.output <- as.data.frame(fect.output)
fect.output$Time <- c(-17:18)
p.fect <- esplot(fect.output,Period = 'Time',Estimate = 'ATT',
                   SE = 'S.E.',CI.lower = "CI.lower", 
                   CI.upper = 'CI.upper',xlim = c(-12,10))
p.fect


## ----hh_impute2, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
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


## ----hh_balance, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
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


## ----load-gs2020, message = FALSE, warning = FALSE----------------------------
data(gs2020)
data <- gs2020
data$cycle <- as.integer(as.numeric(data$cycle/2))
head(data)


## ----gb_panelview_treat, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
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



## ----gb_twfe1, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
model.twfe <- feols(general_sharetotal_A_all ~ cand_A_all + 
                      cand_H_all + cand_B_all | district_final + cycle,
                    data=data, cluster = "district_final") 
summary(model.twfe)


## ----gb_twfe2, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
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


## ----gb_twfeplot, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
twfe.output <- as.data.frame(twfe.est$coeftable[c(1:25),])
twfe.output$Time <- c(c(-16:-2),c(0:9)) + 1 

# plot
p.twfe <- esplot(twfe.output,Period = 'Time',Estimate = 'Estimate',
                               SE = 'Std. Error', xlim = c(-15,1))
p.twfe


## ----gb_pm, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
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



## ----gb_pm1, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
PE.results.pool <- PanelEstimate(PM.results, panel.data = df.pm, pooled = TRUE)
summary(PE.results.pool)


## ----gb_pm2, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE, eval = TRUE----
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


## ----gb_fect, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
model.fect <- fect(Y = "general_sharetotal_A_all", D = "cand_A_all", 
                   X= c("cand_H_all", "cand_B_all"), data = data, 
                   method = "fe", index = index, se = TRUE, 
                   parallel = TRUE, cores = 16, seed = 1234, force = "two-way")

print(model.fect$est.avg)


## ----gb_fectplot, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
fect.output <- as.data.frame(model.fect$est.att)
fect.output$Time <- c(-15:10)
p.fect <- esplot(fect.output,Period = 'Time',Estimate = 'ATT',
                   SE = 'S.E.',CI.lower = "CI.lower", 
                   CI.upper = 'CI.upper', xlim = c(-15,1))
p.fect


## ----gb_fectplot3, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
plot(model.fect)


## ----gb_fectplot4, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
plot(model.fect, type = 'exit')


## ----hh_fectplacebo, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
out.fect.p <- fect(Y = y, X = controls, D = d, data = data, index = index,
                   method = 'fe', se = TRUE, placeboTest = TRUE,
                   placebo.period = c(-2,0))

plot(out.fect.p, proportion = 0.1, stats = "placebo.p")


## ----gb_fectcarryover, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
out.fect.c <- fect(Y = y, X = controls, D = d, data = data, index = index,
                   method = 'fe', se = TRUE, carryoverTest = TRUE, carryover.period = c(1,2))

# plot
plot(out.fect.c,  stats = "carryover.p", ylim = c(-0.15, 0.20))


## ----gb_balance, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
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


## ----wrapper, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
res_st <- did_wrapper(
  data   = hh2019,
  Y      = "nat_rate_ord",
  D      = "indirect",
  index  = c("bfs", "year"),
  method = "st",
  se     = "default"
)
print(res_st)


## ----wrapper_boot, message = FALSE, warning = FALSE, fig.width = 6, fig.height = 4.5, cache=TRUE----
res_st <- did_wrapper(
  data   = hh2019,
  Y      = "nat_rate_ord",
  D      = "indirect",
  index  = c("bfs", "year"),
  method = "st",
  se     = "boot",
  nboots = 200,
  parallel = TRUE
)
print(res_st)


## ----wrapper_plot, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 5, cache=TRUE----
esplot(data = res_st, main = "Stacked DID", xlim = c(-12,10))

