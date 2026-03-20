## Replication file: Liu, Wang, Xu (2022)
# Aim: Generate simulated sample (general treatment pattern)

library("extraDistr")

# e.g. 0 0 1 1 0 1 1 1 -> 0  0  1  2  0  1  2  3
getCumsum <- function(vec){
  b <- unlist(lapply(split(vec,cumsum(c(1, abs(diff(vec))))), cumsum))
  names(b) <- NULL
  return(b)
}

simulateData<-function(
  N,
  TT,
  r, # number of factors
  p, # number of covariates
  beta=c(1,3),
  force=0, # additive fixed effects, 1 = unit, 2 = time, 3 =
  mu=0, # grand mean
  Rtype = "n", ## uniform (u)or normal (n)
  Ftype = "drift", ## factor type: white noise, ar1, drift, trend
  Fsize = 1, ## influence of the factor
  eff.size = 1,
  eff.noise = 1,
  tr.noise = 0.5,
  error.noise = 2,
  factor = NULL,
  seed=NULL
) {

  # panel reshape: T*N -> TN * 1

  if (is.null(seed)==FALSE) {set.seed(seed)}

  #############################
  ## Data generating process
  #############################


  if (r > 0) {
    # loadings
    if (Rtype %in% c("u","uniform")) {
      bound <- sqrt(3)
      lambda <- matrix(runif(N*r, min = 0, max = 2*bound), N, r, byrow=TRUE)
    } else {
      lambda <- matrix(rnorm(N*r, 0.5, 1), N, r, byrow=TRUE)
    }

    # factors
    if (is.null(factor) == TRUE) {
      factor<-matrix(NA,TT,r)
      if (length(Ftype)==1) {Ftype <- rep(Ftype, r)}
      for (i in 1:r) {
        type <- Ftype[i]
        if (type == "ar1") {
          ts <- arima.sim(list(order = c(1,0,0), ar = 0.5), n = (TT+20))[-c(1:20)]
        } else if (type == "drift") {
          ts <- (arima.sim(list(order = c(1,0,0), ar = 0.5), n = (TT+20))[-c(1:20)]
                 + 0.5 * c(1:TT))
        } else if (type == "trend") {
          ts <- c(1:TT) + rnorm(TT, sd = 1)
        } else {
          ts <- rnorm(TT)
        }
        factor[,i] <- ts/sd(ts)
      }
    }
  }




  ## covariates
  if (p>0) {
    X<-array(0,dim=c(TT,N,p))    # regressor matrix, must be T by N by p
    X <- array(rnorm(TT*N*p,sd=1),dim=c(TT,N,p))  # TT*N*p
    beta <- matrix(beta, p, 1)
  }

  # fixed effects
  if (force==1|force==3) { ## unit fixed effect
    alpha <- matrix(rep(rnorm(N,sd=1),each=TT),TT,N)
  }
  if (force==2|force==3) { ## time fixed effect
    drift <- arima.sim(list(order = c(1,0,0), ar = 0.5), n = TT) + 0.5*seq(from=1, length.out = TT, by=0.5)
    xi <- matrix(rep(drift,N),TT,N)
  }

  ## treatment assignment
  D.raw <- matrix(rnorm(TT*N,sd=tr.noise),TT,N)
  if (r > 0){
    D.raw <- D.raw + factor%*%t(lambda) * 5
  }
  if ((force==1|force==3)) {
    D.raw <- D.raw + alpha * 2
  }
  if ((force==2|force==3)) {
    D.raw <- D.raw + xi * 2
  }
  D.raw <- D.raw/10 - 5

  exp.D.raw1 <- exp(D.raw[1,]) # first row
  tr.prob <- matrix(0, TT, N)
  tr.prob[1,] <- exp.D.raw1/(1+exp.D.raw1)
  D <- matrix(0, TT, N)
  D[1,] <- sapply(c(tr.prob[1,]), function(x){rbern(n = 1, prob = x)})
  for (i in 2:TT) {
    exp.D.raw.tmp <- exp(D.raw[i,] + 5 * D[(i-1),])
    tr.prob[i,] <- exp.D.raw.tmp/(1+exp.D.raw.tmp)
    D[i, ] <- sapply(c(tr.prob[i,]), function(x){rbern(n = 1, prob = x)})
  }
  # cumulative
  D.cum <- matrix(0, TT, N)
  for (i in 1:N) {
    D.cum[,i] <- getCumsum(D[,i])
  }



  ## disturbances
  e <- matrix(rnorm(TT*N,sd=error.noise),TT,N)


  ## outcome variable
  Y <- matrix(mu, TT, N)
  if (r>0) {
    if (length(Fsize)==1)  {Fsize <- rep(Fsize,r)}
    for (i in 1:r) {
      factor[,i] <- factor[,i]*Fsize[i]
    }
    Y <- Y + factor%*%t(lambda)
  }
  if (p>0) {
    for (k in 1:p) {
      Y<-Y+X[,,k]*beta[k]
    }
  }

  # fixed effects
  if (force==1|force==3) { ## unit fixed effect
    Y <- Y +3*alpha
  }
  if (force==2|force==3) { ## time fixed effect
    Y <- Y + xi
  }


  ## treatment effect
  eff <- matrix(rnorm(TT*N, mean = eff.size, sd = eff.noise),TT,N) * D.cum
  Y <- Y + eff

  ## disturbance
  Y <- Y + e


  ## panel structure
  panel<-as.data.frame(cbind(
    rep(101:(100+N),each=TT),
    rep(1:TT,N),
    c(Y),
    c(e),
    c(eff),
    c(D.cum),
    c(tr.prob)
    ))
  cname<-c("id","time","Y","error","eff","tr_cum","tr_prob")

  ## treatment indicator
  panel <- cbind(panel,c(D))
  cname <- c(cname,"D")

  ## covar
  if (p>0) {
    for (i in 1:p) {
      panel<-cbind(panel,c(X[,,i]))
      cname<-c(cname,paste("X",i,sep=""))
    }
  }

  ## additive fixed effects
  if (force==1|force==3) {
    panel <- cbind(panel,c(alpha))
    cname <- c(cname,"alpha")
  }
  if (force==2|force==3) {
    panel <- cbind(panel,c(xi))
    cname <- c(cname,"xi")
  }

  if (r>0) {
    for (i in 1:r) {
      panel<-cbind(panel,rep(factor[,i],N))
      cname<-c(cname,paste("F",i,sep=""))
    }
    for (i in 1:r) {
      panel<-cbind(panel,rep(lambda[,i],each=TT))
      cname<-c(cname,paste("L",i,sep=""))
    }
  }
  colnames(panel)<-cname

  if (r>0) for (i in 1:r) {
    panel[,paste("FL",i,sep="")]<-panel[,paste("F",i,sep="")]*panel[,paste("L",i,sep="")]
  }

  return(panel)
}
