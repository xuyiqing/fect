#CIs don't make sense for HonestDiD
library(HonestDiDFEct)
library(fect)

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

out.fect.placebo.ife2 <- fect_sens(
  fect.out     = out.fect.placebo.ife,
  post.periods = 1:10,
  periodMbarvec      = c(1),
)



#Parallel computing does not work
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
                                                                              parallel = TRUE,
                                                                              numPostPeriods = T.post,
                                                                              l_vec = count/sum(count),
                                                                              Mbarvec = seq(0,1,by=0.1))

honest.dte <- HonestDiDFEct::createSensitivityResults_relativeMagnitudes(
  betahat        = beta.hat,
  sigma          = vcov.hat,
  numPrePeriods  = numPrePeriods,
  numPostPeriods = numPostPeriods,
  l_vec          = dte_l,
  Mbarvec        = periodMbarvec,
)
