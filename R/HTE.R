###################################################################
## Heterogeneous Effect Function
###################################################################
#a new function for HTE estimate
fectHTEonce <- function(data,
                        Yname,
                        Dname,
                        Xname = NULL,
                        Moderator,
                        DataType = DataType,
                        Nbins = Nbins,
                        W = NULL, # weight
                        group = NULL, # cohort
                        na.rm = FALSE, # remove missing values
                        index, # c(unit, time) indicators
                        force = "two-way", # fixed effects demeaning
                        r = 0, # number of factors
                        lambda = NULL, # mc method: regularization parameter
                        nlambda = 10, ## mc method: regularization parameter
                        CV = NULL, # cross-validation
                        k = 10, # times of CV
                        cv.prop = 0.1, ## proportion of CV counts
                        cv.treat = FALSE, ## cv targeting treated units
                        cv.nobs = 3,  ## cv taking consecutive units
                        cv.donut = 0, ## cv mspe
                        criterion = "mspe", # for ife model: mspe, pc or both
                        binary = FALSE, # probit model
                        QR = FALSE, # QR or SVD for binary probit
                        method = "fe", # method: e for fixed effects; ife for interactive fe; mc for matrix completion
                        se = FALSE, # report uncertainties
                        vartype = "bootstrap", # bootstrap or jackknife
                        quantile.CI = FALSE,
                        nboots = 200, # number of bootstraps
                        alpha = 0.05, # significance level
                        parallel = TRUE, # parallel computing
                        cores = NULL, # number of cores
                        tol = 0.001, # tolerance level
                        max.iteration = 1000,
                        seed = NULL, # set seed
                        min.T0 = NULL, # minimum T0
                        max.missing = NULL, # maximum missing
                        proportion = 0.3, # use to fit the f test and equivalence test
                        pre.periods = NULL, # fit test period
                        f.threshold = 0.5, # equiv
                        tost.threshold = NULL, # equiv
                        knots = NULL,
                        degree = 2,  # wald = FALSE, # fit test
                        sfe = NULL,
                        cfe = NULL,
                        balance.period = NULL, # the pre and post periods for balanced samples
                        fill.missing = FALSE, # whether to balance missing observations
                        placeboTest = FALSE, # placebo test
                        placebo.period = NULL, # placebo test period
                        carryoverTest = FALSE, # carry-over test
                        carryover.period = NULL, # carry-over period
                        carryover.rm = NULL,
                        loo = FALSE, # leave one period out placebo
                        permute = FALSE, ## permutation test
                        m = 2, ## block length
                        normalize = FALSE, # accelerate option
                        HTE.enp.seq = NULL #parameter used in loess fit estimation
)
{
  out.null<-list(
    method = NULL,
    Y.ct = NULL,
    Y.ct.full = NULL,
    D = NULL,
    Y = NULL,
    X = NULL,
    eff = NULL,
    I = NULL,
    II = NULL,
    att.avg = NA,
    att.avg.boot = NULL,
    att.avg.unit = NA,
    est.avg = c(NA,NA,NA,NA),
    ## supporting
    force = NULL,
    T = NULL,
    N = NULL,
    Ntr = NULL,
    Nco = NULL,
    tr = NULL,
    co = NULL,
    p = NULL,
    r.cv = NULL,
    IC = NULL,
    beta = NULL,
    est = NULL,
    mu = NULL,
    niter = NULL,
    validX = NULL,
    validF = NULL,
    time = NULL,
    att = NULL,
    count = NULL,
    eff.calendar = NULL,
    N.calendar = NULL,
    eff.calendar.fit = NULL,
    calendar.enp = NULL,
    eff.pre = NULL,
    eff.pre.equiv = NULL,
    pre.sd = NULL,
    eff.HTE = NULL,
    Val.HTE = NULL,
    N.HTE = NULL,
    Ntr.HTE = NULL,
    eff.HTE.fit = NULL,
    HTE.enp = NULL,
    bootVal = NULL,
    HTEcoef = NULL,
    KWtest = NULL
    # time.HTE = time.HTE.on,
    # att.HTE = att.HTE.on,
    # count.HTE = count.HTE.on
  )
  HTEvalue = data[,Moderator]

  #split the sample according to the value of moderator variable
  data_list = list()
  if(DataType == "discrete"){    #discrete variable
    HTEuni = unique(as.vector(HTEvalue))
    HTEuni <- sort(HTEuni)
    Val.HTE = HTEuni
    for(i in c(1:length(HTEuni))){
      temp_index = which(data[,Moderator] == HTEuni[i])
      temp_data = data[temp_index,]
      data_list[[i]] = temp_data
    }
  }
  if(DataType == "continuous"){    #continuous variable
    nbins = Nbins
    Val.HTE = rep(NA,nbins)
    quan = 1/nbins
    HTEquantile = quantile(HTEvalue,seq(quan,1,quan))

    for(i in c(1:nbins)){
      if (i == 1){
        temp_index = which((data[,Moderator] < HTEquantile[i]))
      }
      else if (i == nbins){
        temp_index = which((data[,Moderator] >= HTEquantile[i - 1]))
      }
      else {
        temp_index = which((data[,Moderator] >= HTEquantile[i - 1]) & (data[,Moderator] < HTEquantile[i]))
      }
      temp_data = data[temp_index,]
      data_list[[i]] = temp_data
      Val.HTE[i] = quan * i
    }
  }
  out = list()
  kwframe = data.frame()
  id <- index[1]
  time <- index[2]
  for ( i in c(1:length(data_list))){
    temp_data <- data_list[[i]]
    temp_data <- temp_data[order(temp_data[,id],temp_data[,time]),]
    temp_Xname <- Xname
    stable_var<- c()
    if (length(Xname > 0)){
    for (j in 1:length(Xname)){
      if (sum(tapply(temp_data[, Xname[j]], temp_data[,id], var), na.rm = TRUE) == 0){  #deal with the situation that HTE var is unit-invariant in subsample
        stable_var <- c(stable_var, j)
      }
    }
    if (length(stable_var) > 0){
      temp_Xname <- temp_Xname[-stable_var]
    }
    }
    temp_out <- try(fect(data = temp_data,
                         Y = Yname,
                         D = Dname,
                         X = temp_Xname,
                         W = W,
                         group = group,
                         na.rm = na.rm,
                         index = index,
                         force = force,
                         r = r,
                         lambda = lambda,
                         nlambda = nlambda,
                         CV = CV,
                         k = k,
                         cv.prop = cv.prop,
                         cv.treat = cv.treat,
                         cv.nobs = cv.nobs,
                         cv.donut = cv.donut,
                         criterion = criterion,
                         binary = binary,
                         QR = QR,
                         method = method,
                         se = FALSE,    #bootstrap process has been employed in the part above
                         vartype = vartype,
                         quantile.CI = quantile.CI,
                         nboots = nboots,
                         alpha = alpha,
                         parallel = parallel,
                         cores = cores,
                         tol = tol,
                         max.iteration = max.iteration,
                         seed = seed,
                         min.T0 = min.T0,
                         max.missing = max.missing,
                         proportion = proportion,
                         pre.periods = pre.periods,
                         f.threshold = f.threshold,
                         tost.threshold = tost.threshold,
                         knots = knots,
                         degree = degree,
                         sfe = sfe,
                         cfe = cfe,
                         balance.period = balance.period,
                         fill.missing = fill.missing,
                         placeboTest = placeboTest,
                         placebo.period = placebo.period,
                         carryoverTest = carryoverTest,
                         carryover.period = carryover.period,
                         carryover.rm = carryover.rm,
                         loo = loo,
                         permute = permute,
                         m = m,
                         normalize = normalize,
                         Moderator = NULL,
                         HTE.enp.seq = HTE.enp.seq))
    if (class(temp_out) == 'try-error'){
      temp_out = out.null
    }
    Ntr.HTE = sum(temp_data[,Dname], rm.na = TRUE)
    N.HTE = dim(temp_data)[1]
    out[[i]] = temp_out
    out[[i]]$NHTE = N.HTE
    out[[i]]$NtrHTE = Ntr.HTE
    out[[i]]$ValHTE = Val.HTE[i]
  }
  return(out)
}

fectHTEeffonce <- function(eff, HTEvalue, DataType, Nbins, D){
  D.missing <- D
  D.missing[which(D == 0)] <- NA
  out = list()
  if(DataType == "discrete"){   #discrete moderator variable
    HTEuni = unique(as.vector(HTEvalue))
    HTEuni <- sort(HTEuni)
    Val.HTE = HTEuni
    att.avg.HTE <- rep(NA,length(HTEuni))
    N.HTE = rep(0,length(HTEuni))
    Ntr.HTE = rep(0,length(HTEuni))
    for(i in 1:length(HTEuni)){
      INDEX <- D.missing
      INDEX[which(HTEvalue != HTEuni[i])] <- NA
      att.avg.HTE[i] = mean(INDEX * eff, na.rm = TRUE)
      N.HTE[i] = length(which(HTEvalue == HTEuni[i]))
      Ntr.HTE[i] = length(INDEX) - length(is.na(INDEX))
      out[[i]] = list()
      out[[i]]$NHTE = N.HTE[i]
      out[[i]]$NtrHTE = Ntr.HTE[i]
      out[[i]]$att.avg = att.avg.HTE[i]
      out[[i]]$ValHTE = Val.HTE[i]
    }
  }
  else if (DataType == "continuous"){  #continuous variable
    nbins = Nbins
    Val.HTE = rep(NA,nbins)
    quan = 1/nbins
    HTEquantile = quantile(HTEvalue,seq(quan,1,quan))
    att.avg.HTE <- rep(NA,nbins)
    N.HTE = rep(0,nbins)
    Ntr.HTE = rep(0,nbins)
    for(i in c(1:nbins)){
      INDEX <- D.missing
      if (i == 1){
        INDEX[which(HTEvalue >= HTEquantile[i])] <- NA
        N.HTE[i] = length(which(HTEvalue < HTEquantile[i]))
      }
      else if (i == nbins){
        INDEX[which(HTEvalue < HTEquantile[i - 1])] <- NA
        N.HTE[i] = length(which(HTEvalue >= HTEquantile[i - 1]))
      }
      else {
        INDEX[which((HTEvalue >= HTEquantile[i]) | (HTEvalue < HTEquantile[i - 1]))] <- NA
        N.HTE[i] = length(which((HTEvalue >= HTEquantile[i - 1]) & (HTEvalue < HTEquantile[i]) ))
      }
      att.avg.HTE[i] <- mean(INDEX * eff, na.rm = TRUE)
      Val.HTE[i] = quan * i
      Ntr.HTE[i] = sum(INDEX,na.rm = TRUE)
      out[[i]] = list()
      out[[i]]$NHTE = N.HTE[i]
      out[[i]]$NtrHTE = Ntr.HTE[i]
      out[[i]]$att.avg = att.avg.HTE[i]
      out[[i]]$ValHTE = Val.HTE[i]
    }
  }
  return(out)
}
fectHTE <- function(data, # a data frame (long-form)
                    Yname, # outcome
                    Dname, # treatment
                    Xname = NULL, # time-varying covariates
                    W = NULL, # weight
                    group = NULL, # cohort
                    na.rm = FALSE, # remove missing values
                    index, # c(unit, time) indicators
                    force = "two-way", # fixed effects demeaning
                    r = 0, # number of factors
                    lambda = NULL, # mc method: regularization parameter
                    nlambda = 10, ## mc method: regularization parameter
                    CV = NULL, # cross-validation
                    k = 10, # times of CV
                    cv.prop = 0.1, ## proportion of CV counts
                    cv.treat = FALSE, ## cv targeting treated units
                    cv.nobs = 3,  ## cv taking consecutive units
                    cv.donut = 0, ## cv mspe
                    criterion = "mspe", # for ife model: mspe, pc or both
                    binary = FALSE, # probit model
                    QR = FALSE, # QR or SVD for binary probit
                    method = "fe", # method: e for fixed effects; ife for interactive fe; mc for matrix completion
                    se = FALSE, # report uncertainties
                    vartype = "bootstrap", # bootstrap or jackknife
                    quantile.CI = FALSE,
                    nboots = 200, # number of bootstraps
                    alpha = 0.05, # significance level
                    parallel = TRUE, # parallel computing
                    cores = NULL, # number of cores
                    tol = 0.001, # tolerance level
                    max.iteration = 1000,
                    seed = NULL, # set seed
                    min.T0 = NULL, # minimum T0
                    max.missing = NULL, # maximum missing
                    proportion = 0.3, # use to fit the f test and equivalence test
                    pre.periods = NULL, # fit test period
                    f.threshold = 0.5, # equiv
                    tost.threshold = NULL, # equiv
                    knots = NULL,
                    degree = 2,  # wald = FALSE, # fit test
                    sfe = NULL,
                    cfe = NULL,
                    balance.period = NULL, # the pre and post periods for balanced samples
                    fill.missing = FALSE, # whether to balance missing observations
                    placeboTest = FALSE, # placebo test
                    placebo.period = NULL, # placebo test period
                    carryoverTest = FALSE, # carry-over test
                    carryover.period = NULL, # carry-over period
                    carryover.rm = NULL,
                    loo = FALSE, # leave one period out placebo
                    permute = FALSE, ## permutation test
                    m = 2, ## block length
                    normalize = FALSE, # accelerate option
                    Moderator = NULL, #the variable needs heterogeneity estimation
                    DataType = "discrete", #data type of moderator
                    Nbins = NULL, #number of bins
                    HTE.enp.seq = NULL #parameter used in loess fit estimation
) {

  #prepare frequently used variables
  # varnames <- all.vars(formula)
  # Yname <- varnames[1]
  # Dname <- varnames[2]
  out.null<-list(
    method = NULL,
    Y.ct = NULL,
    Y.ct.full = NULL,
    D = NULL,
    Y = NULL,
    X = NULL,
    eff = NULL,
    I = NULL,
    II = NULL,
    att.avg = NA,
    att.avg.boot = NULL,
    att.avg.unit = NA,
    est.avg = c(NA,NA,NA,NA),
    ## supporting
    force = NULL,
    T = NULL,
    N = NULL,
    Ntr = NULL,
    Nco = NULL,
    tr = NULL,
    co = NULL,
    p = NULL,
    r.cv = NULL,
    IC = NULL,
    beta = NULL,
    est = NULL,
    mu = NULL,
    niter = NULL,
    validX = NULL,
    validF = NULL,
    time = NULL,
    att = NULL,
    count = NULL,
    eff.calendar = NULL,
    N.calendar = NULL,
    eff.calendar.fit = NULL,
    calendar.enp = NULL,
    eff.pre = NULL,
    eff.pre.equiv = NULL,
    pre.sd = NULL,
    eff.HTE = NULL,
    Val.HTE = NULL,
    N.HTE = NULL,
    eff.HTE.fit = NULL,
    HTE.enp = NULL,
    bootVal = NULL,
    HTEcoef = NULL,
    KWtest = NULL
    # time.HTE = time.HTE.on,
    # att.HTE = att.HTE.on,
    # count.HTE = count.HTE.on
  )

  # if(DataType == "continuous"){
  #   mode <- 'continuous'
  # }
  # if(DataType == "discrete"){
  #   mode <- 'discrete'
  # }
  # if (length(varnames) > 2) {
  #   Xname <- varnames[3:length(varnames)]
  # } else {
  #   Xname <- NULL
  # }


  if((!Moderator %in% Xname) | (!Moderator %in% colnames(data))){
    stop("Moderator not in X variables or data")
  }

  id = index[1]
  time = index[2]
  data <- data[order(data[,id], data[,time]),]
  TT <- length(unique(data[,time]))
  N <- length(unique(data[,id]))
  D <- matrix(data[, Dname], TT, N)
  Y <- matrix(data[, Yname], TT, N)
  p <- length(Xname)
  X <- array(0, dim = c(TT, N, p))
  id.list <- matrix(data[,id], TT, N)[1,]
  time.list <- matrix(data[,time], TT, N)[,1]
  if (p > 0) {
    for (i in 1:p) {
      X[,,i] <- matrix(data[, Xname[i]], TT, N)
    }
  }
  #mode judgement
  if (DataType == "continuous"){             #continuous moderator
    #uncertainity judgement
    if (se == FALSE){
      result <- fectHTEonce(data = data,Yname = Yname,Dname = Dname,Xname = Xname,Moderator = Moderator, DataType = DataType, Nbins = Nbins, index = index, r = r, force = force, CV = CV)
      return(result)
    }
    else {
      #generate index of bootstrap
      nbins = Nbins
      boot.index <- array(0,dim = c(N,nboots))
      boot.att <- array(NA,dim = c(nbins,nboots))
      sum.D <- colSums(D)
      id.tr <- which(sum.D>0)
      id.co <- which(sum.D==0)
      Nco <- length(id.co)
      Ntr <- length(id.tr)
      I <- matrix(1, TT, N)
      I[is.nan(Y)] <- 0

      out = fectHTEonce(data = data,Yname = Yname,Dname = Dname,Xname = Xname,Moderator = Moderator,DataType = DataType, Nbins = Nbins, index = index, r = r, force = force, CV = CV)

      #initialize the att.boot.list to store the data of att
      att.boot.list = list()
      time.list <- list()
      for(i in 1:length(out)){
        time.list[[i]] = out[[i]]$time
        att.boot.list[[i]] = matrix(NA,length(time.list[[i]]),nboots)
      }

      for(i in 1:nboots){
        repeat {
          fake.co <- sample(id.co,Nco, replace=TRUE)
          if (sum(apply(as.matrix(I[,fake.co]), 1, sum) >= 1) == TT) {
            break
          }
        }
        # fake.co <- sample(id.co,Nco, replace=TRUE)
        repeat {
          fake.tr <- sample(id.tr,Ntr, replace=TRUE)
          if (sum(apply(as.matrix(I[,fake.tr]), 1, sum) >= 1) == TT) {
            break
          }
        }
        # fake.tr <- sample(id.tr,Ntr, replace=TRUE)
        fake.index = c(fake.co,fake.tr)
        boot.index[,i] = fake.index
        temp.index = c()
        temp.id = c()
        for (j in 1:length(fake.index)){
          temp.index.slice = ((fake.index[j] - 1)*TT + 1):(fake.index[j]*TT)
          temp.id.slice = rep(j,TT)
          temp.index = c(temp.index,temp.index.slice)
          temp.id = c(temp.id, temp.id.slice)
        }
        temp.data = data[temp.index,]
        temp.data[[id]] = temp.id
        temp_result = fectHTEonce( data = temp.data,Yname = Yname,Dname = Dname,Xname = Xname,Moderator = Moderator, DataType = DataType, Nbins = Nbins, index = index, r = r, force = force, CV = CV)
        for (k in 1:nbins){
          boot.att[k,i] = temp_result[[k]]$att.avg
          att.boot.list[[k]][,i] = temp_result[[k]]$att[match(time.list[[k]],temp_result[[k]]$time)]
        }
      }

      for (j in 1:nbins){
        temp.boot.rm <- which(is.na(boot.att[j,]))
        temp.avg.att.boot <- boot.att[j,]
        temp.att.boot <- att.boot.list[[j]]
        if(length(temp.boot.rm) > 0){
        temp.avg.att.boot <- boot.att[j,-temp.boot.rm]
        temp.att.boot <- att.boot.list[[j]][,-temp.boot.rm]
        }
        temp.se.avg = sd(temp.avg.att.boot, na.rm = TRUE)
        temp.att.avg = out[[j]]$att.avg
        temp.se.att <- apply(temp.att.boot, 1, function(vec) sd(vec, na.rm=TRUE))
        temp.att = out[[j]]$att
        if(quantile.CI == FALSE){
          temp.CI.avg <- c(temp.att.avg - temp.se.avg * qnorm(1-alpha/2), temp.att.avg + temp.se.avg * qnorm(1-alpha/2))
          temp.pvalue.avg <- (1-pnorm(abs(temp.att.avg/temp.se.avg)))*2

          temp.CI.att <- cbind(temp.att - temp.se.att * qnorm(1-alpha/2), temp.att + temp.att.avg * qnorm(1-alpha/2))
          temp.pvalue.att <- (1-pnorm(abs(temp.att/temp.se.att)))*2
        }
        else{
          temp.CI.avg <- quantile(temp.avg.att.boot,c(alpha/2,1-alpha/2), na.rm=TRUE)
          pvalue.avg  <- get.pvalue(temp.avg.att.boot)

          temp.CI.att <- t(apply(temp.att.boot, 1, function(vec) quantile(vec,c(alpha/2, 1 - alpha/2), na.rm=TRUE)))
          pvalue.att   <- apply(temp.att.boot, 1, get.pvalue)
        }
        temp.est.avg <- t(as.matrix(c(temp.att.avg, temp.se.avg, temp.CI.avg, temp.pvalue.avg)))
        colnames(temp.est.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")

        temp.est.att <- cbind(temp.att, temp.se.att, temp.CI.att, temp.pvalue.att, out[[j]]$count)
        colnames(temp.est.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper","p.value", "count")
        rownames(temp.est.att) <- out[[j]]$time

        out[[j]]$est.avg = temp.est.avg
        out[[j]]$est.att = temp.est.att
      }
      return(out)
    }
  }
  else if (DataType == "discrete"){             #discrete moderator
    #uncertainity judgement
    if (se == FALSE){
      result <- fectHTEonce(data = data,Yname = Yname, Dname= Dname,Xname = Xname,Moderator = Moderator, DataType = DataType, Nbins = Nbins, index = index, r = r, force = force, CV = CV)
      return(result)
    }
    else {
      #generate index of bootstrap
      HTEvalue = data[,Moderator]
      HTEuni = unique(as.vector(HTEvalue))
      nbins = length(HTEuni)
      boot.index <- array(0,dim = c(N,nboots))
      boot.att <- array(NA,dim = c(nbins,nboots))
      sum.D <- colSums(D)
      id.tr <- which(sum.D>0)
      id.co <- which(sum.D==0)
      Nco <- length(id.co)
      Ntr <- length(id.tr)
      I <- matrix(1, TT, N)
      I[is.nan(Y)] <- 0

      out = fectHTEonce(data = data,Yname = Yname, Dname = Dname, Xname = Xname,Moderator = Moderator,DataType = DataType, Nbins = Nbins, index = index, r = r, force = force, CV = CV)

      #initialize the att.boot.list to store the data of att
      att.boot.list = list()
      time.list <- list()
      for(i in 1:length(out)){
        time.list[[i]] = out[[i]]$time
        att.boot.list[[i]] = matrix(NA,length(time.list[[i]]),nboots)
      }

      for(i in 1:nboots){
        repeat {
          fake.co <- sample(id.co, Nco, replace=TRUE)
          if (sum(apply(as.matrix(I[,fake.co]), 1, sum) >= 1) == TT) {
            break
          }
        }
        # fake.co <- sample(id.co,Nco, replace=TRUE)
        repeat {
          fake.tr <- sample(id.tr,Ntr, replace=TRUE)
          if (sum(apply(as.matrix(I[,fake.tr]), 1, sum) >= 1) == TT) {
            break
          }
        }
        # fake.tr <- sample(id.tr,Ntr, replace=TRUE)
        fake.index = c(fake.co,fake.tr)
        boot.index[,i] = fake.index
        temp.index = c()
        temp.id = c()
        for (j in 1:length(fake.index)){
          temp.index.slice = ((fake.index[j] - 1)*TT + 1):(fake.index[j]*TT)
          temp.id.slice = rep(j,TT)
          temp.index = c(temp.index,temp.index.slice)
          temp.id = c(temp.id, temp.id.slice)
        }
        temp.data = data[temp.index,]
        temp.data[[id]] = temp.id
        temp_result = fectHTEonce( data = temp.data,Yname = Yname,Dname = Dname,Xname = Xname,Moderator = Moderator, DataType = DataType, Nbins = Nbins, index = index, r = r, force = force, CV = CV)
        for (k in 1:nbins){
          boot.att[k,i] = temp_result[[k]]$att.avg
          att.boot.list[[k]][,i] = temp_result[[k]]$att[match(time.list[[k]],temp_result[[k]]$time)]
        }
      }

      for (j in 1:nbins){
        temp.boot.rm <- which(is.na(boot.att[j,]))
        temp.avg.att.boot <- boot.att[j,]
        temp.att.boot <- att.boot.list[[j]]
        if(length(temp.boot.rm) > 0){
        temp.avg.att.boot <- boot.att[j,-temp.boot.rm]
        temp.att.boot <- att.boot.list[[j]][,-temp.boot.rm]
        }
        temp.se.avg = sd(temp.avg.att.boot, na.rm = TRUE)
        temp.att.avg = out[[j]]$att.avg
        temp.se.att <- apply(temp.att.boot, 1, function(vec) sd(vec, na.rm=TRUE))
        temp.att = out[[j]]$att
        if(quantile.CI == FALSE){
          temp.CI.avg <- c(temp.att.avg - temp.se.avg * qnorm(1-alpha/2), temp.att.avg + temp.se.avg * qnorm(1-alpha/2))
          temp.pvalue.avg <- (1-pnorm(abs(temp.att.avg/temp.se.avg)))*2

          temp.CI.att <- cbind(temp.att - temp.se.att * qnorm(1-alpha/2), temp.att + temp.att.avg * qnorm(1-alpha/2))
          temp.pvalue.att <- (1-pnorm(abs(temp.att/temp.se.att)))*2
        }
        else{
          temp.CI.avg <- quantile(temp.avg.att.boot,c(alpha/2,1-alpha/2), na.rm=TRUE)
          pvalue.avg  <- get.pvalue(temp.avg.att.boot)

          temp.CI.att <- t(apply(temp.att.boot, 1, function(vec) quantile(vec,c(alpha/2, 1 - alpha/2), na.rm=TRUE)))
          pvalue.att   <- apply(temp.att.boot, 1, get.pvalue)
        }
        temp.est.avg <- t(as.matrix(c(temp.att.avg, temp.se.avg, temp.CI.avg, temp.pvalue.avg)))
        colnames(temp.est.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")

        temp.est.att <- cbind(temp.att, temp.se.att, temp.CI.att, temp.pvalue.att, out[[j]]$count)
        colnames(temp.est.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper","p.value", "count")
        rownames(temp.est.att) <- out[[j]]$time

        out[[j]]$est.avg = temp.est.avg
        out[[j]]$est.att = temp.est.att
      }
      return(out)
    }
  }


}


testHTE <- function(out){
  kwframe = data.frame()
  for(i in 1:length(out)){
    if(! is.null(out[[i]]$att.avg.boot)){
      temp_dataframe = data.frame(val = c(out[[i]]$att.avg.boot),group = rep(i,length(out[[i]]$att.avg.boot)))
      kwframe = rbind(kwframe,temp_dataframe)
    }
  }
  kwtest = NULL
  aovtest = NULL
  if (! is.null(kwframe)){
    kwtest = list(kruskal.test(val~group,data = kwframe))
    aovtest = list(aov(val~group,data = kwframe))
  }
  result = list(kwtest = kwtest,aovtest = aovtest)
  return(result)
}






plotHTE <- function(x,
                      version = "whole", # whole, seperate
                      mode = "discrete", # discrete, continuous
                      type = 'avg', # avg, dynamic
                      loo = FALSE,
                      highlight = NULL, ## for carryover test and placebo test
                      plot.ci = "0.95", ## "0.9", "0.95", "none"
                      show.points = NULL,
                      show.group = NULL,
                      bound = NULL, # "none", "min", "equiv", "both"
                      vis = NULL,
                      count = TRUE,
                      proportion = 0.3, # control the xlim
                      pre.periods = NULL, # for testing
                      f.threshold = NULL, # equiv f
                      tost.threshold = NULL, # pre-trend placebo carryover
                      effect.bound.ratio = FALSE,
                      stats = NULL,       ## "none", "F.p", "F.equiv.p", "placebo.p", "carryover.p", "equiv.p"
                      stats.labs = NULL,
                      raw = "none",    ## "none", "band", "all"
                      main = NULL,
                      xlim = NULL,
                      ylim = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      gridOff = FALSE,
                      legendOff = FALSE,
                      legend.pos = NULL,
                      legend.nrow = NULL,
                      legend.labs = NULL,
                      stats.pos = NULL,
                      theme.bw = TRUE,
                      nfactors = NULL,
                      include.FE = TRUE,
                      id = NULL,
                      cex.main = NULL,
                      cex.main.sub = NULL,
                      cex.axis = NULL,
                      cex.lab = NULL,
                      cex.legend = NULL,
                      cex.text = NULL,
                      axis.adjust = FALSE,
                      axis.lab = "both",
                      axis.lab.gap = c(0, 0),
                      shade.post = FALSE,
                      start0 = FALSE,
                      return.test = FALSE,
                      balance = NULL,
                      weight = NULL,
                      save_path = NULL
                      ){
  # stats <- "none"
  # # names for all statistics
  # if (!("none" %in% stats)) {
  #   if (is.null(stats.labs)==FALSE) {
  #     if (length(stats.labs)!=length(stats)) {
  #       stop("\"stats.lab\" should have the same length as \"stats\".")
  #     }
  #   }
  #   else {
  #     stats.labs <- rep(NA, length(stats))
  #     for (i in 1:length(stats)) {
  #       if (stats[i] == "F.p") {
  #         stats.labs[i] <- "F test p-value"
  #       }
  #       if (stats[i] == "F.equiv.p") {
  #         stats.labs[i] <- "F equivalence test p-value"
  #       }
  #       if (stats[i] == "F.stat") {
  #         stats.labs[i] <- "F statistics"
  #       }
  #       if (stats[i] == "placebo.p") {
  #         stats.labs[i] <- "Placebo test p-value"
  #       }
  #       if (stats[i] == "carryover.p") {
  #         stats.labs[i] <- "Carryover effect test p-value"
  #       }
  #       if (stats[i] == "equiv.p") {
  #         if(placeboTest){
  #           stats.labs[i] <- "Placebo equivalence test p-value"
  #         }
  #         else if(carryoverTest){
  #           stats.labs[i] <- "Carryover effect equivalence test p-value"
  #         }
  #         else{
  #           stats.labs[i] <- "Equivalence test p-value"
  #         }
  #       }
  #     }
  #   }
  # }
  #
  # titles; xlim and ylim
  ytitle <- NULL
  bound.old <- bound
  maintext <- "ATT by Moderator"
  ytitle <- paste("Effect on Dependent Variable")
  #### font size
  ## title
  if (is.null(cex.main)==FALSE) {
    if (is.numeric(cex.main)==FALSE) {
      stop("\"cex.main\" is not numeric.")
    }
    cex.main <- 16 * cex.main
  } else {
    cex.main <- 16
  }
  ## subtitle
  if (is.null(cex.main.sub)==FALSE) {
    if (is.numeric(cex.main.sub)==FALSE) {
      stop("\"cex.main.sub\" is not numeric.")
    }
    cex.main.sub <- 16 * cex.main.sub
  } else {
    cex.main.sub <- 16
  }
  ## axis label
  if (is.null(cex.lab)==FALSE) {
    if (is.numeric(cex.lab)==FALSE) {
      stop("\"cex.lab\" is not numeric.")
    }
    cex.lab <- 15 * cex.lab
  } else {
    cex.lab <- 15
  }
  ## axis number
  if (is.null(cex.axis)==FALSE) {
    if (is.numeric(cex.axis)==FALSE) {
      stop("\"cex.axis\" is not numeric.")
    }
    cex.axis <- 15 * cex.axis
  }  else {
    cex.axis <- 15
  }
  ## legend
  if (is.null(cex.legend)==FALSE) {
    if (is.numeric(cex.legend)==FALSE) {
      stop("\"cex.legend\" is not numeric.")
    }
    cex.legend <- 15 * cex.legend
  }  else {
    cex.legend <- 15
  }
  ## text
  if (is.null(cex.text)==FALSE) {
    if (is.numeric(cex.text)==FALSE) {
      stop("\"cex.text\" is not numeric.")
    }
    cex.text <- 5 * cex.text
  }  else {
    cex.text <- 5
  }

  ## text label position
  if (!is.null(stats.pos)) {
    if (length(stats.pos) != 2) {
      stop(" \"stats.pos\" must be of length 2. ")
    }
  }
  if(type == "avg" & version == 'seperate'){
  CI <- NULL  #decide if we need plot uncertainties
  for(i in 1:length(x)){
    if ("est.avg" %in% names(x[[i]])){
      if(sum(is.na(x[[i]]$est.avg)) != 4){
        CI <- TRUE
      }
    }
  }
  if (is.null(CI)) {
    CI <- FALSE
  }
  if(plot.ci=="none"){
    CI <- FALSE
  }
  ## axes labels
  if (is.null(xlab) == TRUE) {
    xlab <- "Moderator"
  } else if (xlab == "") {
    xlab <- NULL
  }

  if (is.null(ylab) == TRUE) {
    ylab <- ytitle
  } else if (ylab == "") {
    ylab <- NULL
  }

  ## y=0 line type
  lcolor <- "white"
  lwidth <- 2
  if (theme.bw == TRUE) {
    lcolor <- "#AAAAAA70"
    lwidth <- 1.5
  }

  if (CI == FALSE) {
    message("Uncertainty estimates not available.\n")
    YVAR = rep(NA,length(x))
    XVAR = rep(NA,length(x))
    for (i in 1:length(x)){
      YVAR[i] = x[[i]]$att.avg
      XVAR[i] = x[[i]]$ValHTE
    }
    data.toplot.main <- cbind.data.frame(YVAR,XVAR)

    if (length(ylim) != 0) {
      rect.length <- (ylim[2] - ylim[1]) / 5
      rect.min <- ylim[1]
    }
    else {
      rect.length <- (max(c(YVAR), na.rm = TRUE) - min(c(YVAR), na.rm = TRUE))/2
      rect.min <- min(c(YVAR), na.rm = TRUE) - rect.length
    }
  }
  else {
  YVAR = rep(NA,length(x))
  XVAR = rep(NA,length(x))
  YMIN = rep(NA,length(x))
  YMAX = rep(NA,length(x))
  for (i in 1:length(x)){
    YVAR[i] = x[[i]]$est.avg[1]
    YMIN[i] = x[[i]]$est.avg[3]
    YMAX[i] = x[[i]]$est.avg[4]
    XVAR[i] = x[[i]]$ValHTE
  }
  data.toplot.main <- cbind.data.frame(YVAR,YMIN,YMAX,XVAR)
  if (length(ylim) != 0) {
    rect.length <- (ylim[2] - ylim[1]) / 5
    rect.min <- ylim[1]
  } else {
    rect.length <- (max(c(YMAX), na.rm = TRUE) - min(c(YMIN), na.rm = TRUE))/2
    rect.min <- min(c(YMIN), na.rm = TRUE) - rect.length
  }
 }

p <- ggplot(data.toplot.main)

## xlab and ylab
p <- p + xlab(xlab) +  ylab(ylab)
if(x[[1]]$DataType == "continuous"){
  x_labels = rep("",length(XVAR))
  for(i in 1:length(XVAR)){
    if (i == 1){
      x_labels[i] = paste("0%-",round(100*XVAR[i],1),"%")
    }
    else{
      x_labels[i] = paste(round(100*XVAR[i-1],1),"%-",round(100*XVAR[i],1),"%")
    }

  }
  p <- p + scale_x_continuous(breaks = XVAR, labels = x_labels)
}
## theme
if (theme.bw == TRUE) {
  p <- p + theme_bw()
}

## grid
if (gridOff == TRUE) {
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

# horizontal 0 line
#p <- p + geom_hline(yintercept = 0, colour = lcolor,size = lwidth)


if(CI==FALSE){
  #p <- p + geom_point(aes(x=XVAR,y=YVAR,color='gray50',fill='gray50',alpha=1,size=1.2))
  p <- p + geom_point(aes(x=XVAR,y=YVAR),color='gray50',fill='gray50',alpha=1)
} else {
  #p <- p + geom_pointrange(aes(x=XVAR,y=YVAR,ymin=YMIN,ymax=YMAX,color='gray50',fill='gray50',alpha=1,size=0.6))
  p <- p + geom_pointrange(aes(x=XVAR,y=YVAR,ymin=YMIN,ymax=YMAX),color='gray50',fill='gray50',alpha=1)
}
if(count==TRUE){
  NcoCOUNT = rep(NA,length(x))
  NtrCOUNT = rep(NA,length(x))
  NCOUNT = rep(NA,length(x))
  for (i in 1:length(x)){
    NtrCOUNT[i] = x[[i]]$NtrHTE
    NcoCOUNT[i] = x[[i]]$NHTE - NtrCOUNT[i]
    NCOUNT[i] = x[[i]]$NHTE
  }
  T.start <- c()
  T.end <- c()
  ymin <- c()
  ymaxco <- c()
  ymaxtr <- c()
  T.gap <- (max(XVAR)-min(XVAR))/length(XVAR)
  for(i in c(1:length(XVAR))){
    T.start <- c(T.start,XVAR[i]-0.25*T.gap)
    T.end <- c(T.end,XVAR[i]+0.25*T.gap)
    ymin <- c(ymin, rect.min)
    ymaxco <- c(ymaxco, rect.min+rect.length*NcoCOUNT[i]/max(NCOUNT))
    ymaxtr <- c(ymaxtr, rect.min+rect.length*NCOUNT[i]/max(NCOUNT))
  }
  data.toplotco <- cbind.data.frame(xmin=T.start,
                                  xmax=T.end,
                                  ymin=ymin,
                                  ymax=ymaxco)
  data.toplottr <- cbind.data.frame(xmin=T.start,
                                    xmax=T.end,
                                    ymin=ymaxco,
                                    ymax=ymaxtr)
  max.count.pos <- mean(XVAR[which.max(NCOUNT)])
  p <- p + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),data=data.toplotco,fill='gray50',alpha=0.3,size=0.3,color='black')
  p <- p + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),data=data.toplottr,fill='red',alpha=0.3,size=0.3,color='black')
  p <- p + annotate("text", x = max.count.pos - 0.02 * T.gap,
                    y = max(data.toplottr$ymax) + 0.2 * rect.length,
                    label = max(NCOUNT), size = cex.text * 0.8, hjust = 0.5)
}
if (is.null(main) == TRUE) {
  p <- p + ggtitle(maintext) + theme(plot.title = element_text(hjust = 0.5))
} else if (main!=""){
  p <- p + ggtitle(main)
}

if (is.null(save_path) == FALSE){
  filename <- file.path(save_path, paste0("plot_avg.png"))
  ggsave(filename,p)
  }

  }
  else if(type == 'dynamic' & version == 'seperate'){
    CI <- NULL  #decide if we need plot uncertainties
    if ("est.att" %in% names(x[[1]])){
      if(! is.null(x[[1]]$est.att)){
        CI <- TRUE
      }
    }
    if (is.null(CI)) {
      CI <- FALSE
    }
    if(plot.ci=="none"){
      CI <- FALSE
    }
    p_collection = list()

    for (i in 1:length(x)){
      ## axes labels
      if (is.null(xlab) == TRUE) {
        xlab <- "Time since the Treatment began"
      }
      else if (xlab == "") {
        xlab <- NULL
      }

      if (is.null(ylab) == TRUE) {
        ylab <- ytitle
      }
      else if (ylab == "") {
        ylab <- NULL
      }

      ## y=0 line type
      lcolor <- "white"
      lwidth <- 2
      if (theme.bw == TRUE) {
        lcolor <- "#AAAAAA70"
        lwidth <- 1.5
      }

      if (CI == FALSE) {
        message("Uncertainty estimates not available.\n")
        temp.att <- x[[i]]$att
        temp.time <- x[[i]]$time
        temp.count <- x[[i]]$count
        YVAR = rep(NA,length(temp.att))
        XVAR = rep(NA,length(temp.time))
        NCOUNT = rep(NA,length(temp.count))
        for (j in 1:length(temp.att)){
          YVAR[j] = temp.att[j]
          XVAR[j] = temp.time[j]
          NCOUNT[j] = temp.count[j]
        }
        data.toplot.main <- cbind.data.frame(YVAR,XVAR)

        if (length(ylim) != 0) {
          rect.length <- (ylim[2] - ylim[1]) / 5
          rect.min <- ylim[1]
        } else {
          rect.length <- (max(c(YVAR), na.rm = TRUE) - min(c(YVAR), na.rm = TRUE))/2
          rect.min <- min(c(YVAR), na.rm = TRUE) - rect.length
        }
      } else {
        temp.est.att = x[[i]]$est.att
        temp.time <- x[[i]]$time
        rm.pos <- which(is.na(temp.est.att[,3]))
        if (length(rm.pos) != 0){
          temp.est.att <- temp.est.att[-rm.pos,]
          temp.time <- x$time.HTE[[i]][-rm.pos]
        }
        YVAR = rep(NA,dim(temp.est.att)[1])
        XVAR = rep(NA,dim(temp.est.att)[1])
        YMIN = rep(NA,dim(temp.est.att)[1])
        YMAX = rep(NA,dim(temp.est.att)[1])
        NCOUNT = rep(NA,dim(temp.est.att)[1])
        for (j in 1:dim(temp.est.att)[1]){
          YVAR[j] = temp.est.att[j,1]
          YMIN[j] = temp.est.att[j,3]
          YMAX[j] = temp.est.att[j,4]
          XVAR[j] = temp.time[j]
          NCOUNT[j] = temp.est.att[j,6]
        }
        data.toplot.main <- cbind.data.frame(YVAR,YMIN,YMAX,XVAR)
        if (length(ylim) != 0) {
          rect.length <- (ylim[2] - ylim[1]) / 5
          rect.min <- ylim[1]
        } else {
          rect.length <- (max(c(YMAX), na.rm = TRUE) - min(c(YMIN), na.rm = TRUE))/2
          rect.min <- min(c(YMIN), na.rm = TRUE) - rect.length
        }
      }

      temp.p <- ggplot(data.toplot.main)
      ## xlab and ylab
      temp.p <- temp.p + xlab(xlab) +  ylab(ylab)

      ## theme
      if (theme.bw == TRUE) {
        temp.p <- temp.p + theme_bw()
      }

      ## grid
      if (gridOff == TRUE) {
        temp.p <- temp.p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      }

      # horizontal 0 line
      #temp.p <- temp.p + geom_hline(yintercetemp.pt = 0, colour = lcolor,size = lwidth)


      if(CI==FALSE){
        #temp.p <- temp.p + geom_point(aes(x=XVAR,y=YVAR,color='gray50',fill='gray50',alpha=1,size=1.2))
        temp.p <- temp.p + geom_point(aes(x=XVAR,y=YVAR),color='gray50',fill='gray50',alpha=1)
      } else {
        #temp.p <- temp.p + geom_pointrange(aes(x=XVAR,y=YVAR,ymin=YMIN,ymax=YMAX,color='gray50',fill='gray50',alpha=1,size=0.6))
        temp.p <- temp.p + geom_pointrange(aes(x=XVAR,y=YVAR,ymin=YMIN,ymax=YMAX),color='gray50',fill='gray50',alpha=1)
      }
      if(count==TRUE){

        T.start <- c()
        T.end <- c()
        ymin <- c()
        ymax <- c()
        T.gap <- (max(XVAR)-min(XVAR))/length(XVAR)
        for(j in c(1:length(XVAR))){
          T.start <- c(T.start,XVAR[j]-0.25*T.gap)
          T.end <- c(T.end,XVAR[j]+0.25*T.gap)
          ymin <- c(ymin, rect.min)
          ymax <- c(ymax, rect.min+rect.length*NCOUNT[j]/max(NCOUNT))
        }
        data.toplot <- cbind.data.frame(xmin=T.start,
                                        xmax=T.end,
                                        ymin=ymin,
                                        ymax=ymax)
        max.count.pos <- mean(XVAR[which.max(NCOUNT)])
        temp.p <- temp.p + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),data=data.toplot,fill='gray50',alpha=0.3,size=0.3,color='black')
        temp.p <- temp.p + annotate("text", x = max.count.pos - 0.02 * T.gap,
                                    y = max(data.toplot$ymax) + 0.2 * rect.length,
                                    label = max(NCOUNT), size = cex.text * 0.8, hjust = 0.5)
      }
      ## title
      if (is.null(main) == TRUE) {
        temp.title = paste0("Dynamic Effect when Moderator = ",x[[i]]$ValHTE)
        temp.p <- temp.p + ggtitle(temp.title) + theme(plot.title = element_text(hjust = 0.5))
      } else if (main!=""){
        temp.p <- temp.p + ggtitle(main)
      }

      p_collection[[i]] = temp.p

      if (is.null(save_path) == FALSE){
        filename <- file.path(save_path, paste0("plot_dynamic_",i,".png"))
        ggsave(filename,temp.p)
      }

    }
    p <- plot_grid(plotlist = p_collection, ncol = 1)
  }

  else if(type == "avg" & version == 'whole'){
    CI <- NULL  #decide if we need plot uncertainties
    if ("est.avg.HTE" %in% names(x)){
      CI <- TRUE
    }
    if (is.null(CI)) {
      CI <- FALSE
    }
    if(plot.ci=="none"){
      CI <- FALSE
    }
    ## axes labels
    if (is.null(xlab) == TRUE) {
      xlab <- "Moderator"
    }
    else if (xlab == "") {
      xlab <- NULL
    }

    if (is.null(ylab) == TRUE) {
      ylab <- ytitle
    }
    else if (ylab == "") {
      ylab <- NULL
    }

    ## y=0 line type
    lcolor <- "white"
    lwidth <- 2
    if (theme.bw == TRUE) {
      lcolor <- "#AAAAAA70"
      lwidth <- 1.5
    }

    if (CI == FALSE) {
      message("Uncertainty estimates not available.\n")
      YVAR = rep(NA,length(x$avg.HTE))
      XVAR = rep(NA,length(x$avg.HTE))
      for (i in 1:length(x$avg.HTE)){
        YVAR[i] = x$avg.HTE[i]
        XVAR[i] = x$Val.HTE[i]
      }
      data.toplot.main <- cbind.data.frame(YVAR,XVAR)

      if (length(ylim) != 0) {
        rect.length <- (ylim[2] - ylim[1]) / 5
        rect.min <- ylim[1]
      } else {
        rect.length <- (max(c(YVAR), na.rm = TRUE) - min(c(YVAR), na.rm = TRUE))/2
        rect.min <- min(c(YVAR), na.rm = TRUE) - rect.length
      }
    }
    else {
      YVAR = rep(NA,length(x$avg.HTE))
      XVAR = rep(NA,length(x$avg.HTE))
      YMIN = rep(NA,length(x$avg.HTE))
      YMAX = rep(NA,length(x$avg.HTE))
      for (i in 1:length(x$avg.HTE)){
        YVAR[i] = x$est.avg.HTE[i,1]
        YMIN[i] = x$est.avg.HTE[i,3]
        YMAX[i] = x$est.avg.HTE[i,4]
        XVAR[i] = x$Val.HTE[i]
      }
      data.toplot.main <- cbind.data.frame(YVAR,YMIN,YMAX,XVAR)
      if (length(ylim) != 0) {
        rect.length <- (ylim[2] - ylim[1]) / 5
        rect.min <- ylim[1]
      } else {
        rect.length <- (max(c(YMAX), na.rm = TRUE) - min(c(YMIN), na.rm = TRUE))/2
        rect.min <- min(c(YMIN), na.rm = TRUE) - rect.length
      }
    }

    p <- ggplot(data.toplot.main)
    ## xlab and ylab
    p <- p + xlab(xlab) +  ylab(ylab)
    if(x$DataType == "continuous"){
      x_labels = rep("",length(XVAR))
      for(i in 1:length(XVAR)){
        if (i == 1){
          x_labels[i] = paste("0%-",round(100*XVAR[i],1),"%")
        }
        else{
          x_labels[i] = paste(round(100*XVAR[i-1],1),"%-",round(100*XVAR[i],1),"%")
        }

      }
      p <- p + scale_x_continuous(breaks = XVAR, labels = x_labels)
    }
    ## theme
    if (theme.bw == TRUE) {
      p <- p + theme_bw()
    }

    ## grid
    if (gridOff == TRUE) {
      p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }

    # horizontal 0 line
    #p <- p + geom_hline(yintercept = 0, colour = lcolor,size = lwidth)


    if(CI==FALSE){
      #p <- p + geom_point(aes(x=XVAR,y=YVAR,color='gray50',fill='gray50',alpha=1,size=1.2))
      p <- p + geom_point(aes(x=XVAR,y=YVAR),color='gray50',fill='gray50',alpha=1)
    } else {
      #p <- p + geom_pointrange(aes(x=XVAR,y=YVAR,ymin=YMIN,ymax=YMAX,color='gray50',fill='gray50',alpha=1,size=0.6))
      p <- p + geom_pointrange(aes(x=XVAR,y=YVAR,ymin=YMIN,ymax=YMAX),color='gray50',fill='gray50',alpha=1)
    }
    if(count==TRUE){
      NcoCOUNT = rep(NA,length(x$avg.HTE))
      NtrCOUNT = rep(NA,length(x$avg.HTE))
      NCOUNT = rep(NA,length(x$avg.HTE))
      for (i in 1:length(x$avg.HTE)){
        #    NtrCOUNT[i] = x[[i]]$NtrHTE
        #   NcoCOUNT[i] = x[[i]]$NHTE - NtrCOUNT[i]
        NCOUNT[i] = x$N.HTE[i]
        NtrCOUNT[i] = x$Ntr.HTE[i]
        NcoCOUNT[i] = NCOUNT[i] - NtrCOUNT[i]
      }
      T.start <- c()
      T.end <- c()
      ymin <- c()
      ymaxco <- c()
      ymaxtr <- c()
      T.gap <- (max(XVAR)-min(XVAR))/length(XVAR)
      for(i in c(1:length(XVAR))){
        T.start <- c(T.start,XVAR[i]-0.25*T.gap)
        T.end <- c(T.end,XVAR[i]+0.25*T.gap)
        ymin <- c(ymin, rect.min)
        ymaxco <- c(ymaxco, rect.min+rect.length*NcoCOUNT[i]/max(NCOUNT))
        ymaxtr <- c(ymaxtr, rect.min+rect.length*NCOUNT[i]/max(NCOUNT))
      }
      data.toplotco <- cbind.data.frame(xmin=T.start,
                                        xmax=T.end,
                                        ymin=ymin,
                                        ymax=ymaxco)
      data.toplottr <- cbind.data.frame(xmin=T.start,
                                        xmax=T.end,
                                        ymin=ymaxco,
                                        ymax=ymaxtr)
      max.count.pos <- mean(XVAR[which.max(NCOUNT)])
      p <- p + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),data=data.toplotco,fill='gray50',alpha=0.3,size=0.3,color='black')
      p <- p + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),data=data.toplottr,fill='red',alpha=0.3,size=0.3,color='black')
      p <- p + annotate("text", x = max.count.pos - 0.02 * T.gap,
                        y = max(data.toplottr$ymax) + 0.2 * rect.length,
                        label = max(NCOUNT), size = cex.text * 0.8, hjust = 0.5)
    }
    if (is.null(main) == TRUE) {
      p <- p + ggtitle(maintext) + theme(plot.title = element_text(hjust = 0.5))
    } else if (main!=""){
      p <- p + ggtitle(main)
    }

    if (is.null(save_path) == FALSE){
      filename <- file.path(save_path, paste0("plot_avg.png"))
      ggsave(filename,p)
    }
  }

  else if(type == 'dynamic' & version == 'whole'){
    CI <- NULL  #decide if we need plot uncertainties
    if ("est.att.HTE" %in% names(x)){
      if(! is.null(x$est.att.HTE)){
        CI <- TRUE
      }
    }
    if (is.null(CI)) {
      CI <- FALSE
    }
    if(plot.ci=="none"){
      CI <- FALSE
    }
    p_collection = list()

      for (i in 1:length(x$att.HTE)){
        ## axes labels
        if (is.null(xlab) == TRUE) {
          xlab <- "Moderator"
        }
        else if (xlab == "") {
          xlab <- NULL
        }

        if (is.null(ylab) == TRUE) {
          ylab <- ytitle
        }
        else if (ylab == "") {
          ylab <- NULL
        }

        ## y=0 line type
        lcolor <- "white"
        lwidth <- 2
        if (theme.bw == TRUE) {
          lcolor <- "#AAAAAA70"
          lwidth <- 1.5
        }

        if (CI == FALSE) {
          message("Uncertainty estimates not available.\n")
          temp.att <- x$att.HTE[[i]]
          temp.time <- x$time.HTE[[i]]
          temp.count <- x$count.HTE[[i]]
          YVAR = rep(NA,length(temp.att))
          XVAR = rep(NA,length(temp.time))
          NCOUNT = rep(NA,length(temp.count))
          for (j in 1:length(temp.att)){
            YVAR[j] = temp.att[j]
            XVAR[j] = temp.time[j]
            NCOUNT[j] = temp.count[j]
          }
          data.toplot.main <- cbind.data.frame(YVAR,XVAR)

          if (length(ylim) != 0) {
            rect.length <- (ylim[2] - ylim[1]) / 5
            rect.min <- ylim[1]
          } else {
            rect.length <- (max(c(YVAR), na.rm = TRUE) - min(c(YVAR), na.rm = TRUE))/2
            rect.min <- min(c(YVAR), na.rm = TRUE) - rect.length
          }
        }
        else {
          temp.est.att = x$est.att.HTE[[i]]
          temp.time <- x$time.HTE[[i]]
          rm.pos <- which(is.na(temp.est.att[,3]))
          if (length(rm.pos) != 0){
            temp.est.att <- temp.est.att[-rm.pos,]
            temp.time <- x$time.HTE[[i]][-rm.pos]
            }
          YVAR = rep(NA,dim(temp.est.att)[1])
          XVAR = rep(NA,dim(temp.est.att)[1])
          YMIN = rep(NA,dim(temp.est.att)[1])
          YMAX = rep(NA,dim(temp.est.att)[1])
          NCOUNT = rep(NA,dim(temp.est.att)[1])
          for (j in 1:dim(temp.est.att)[1]){
            YVAR[j] = temp.est.att[j,1]
            YMIN[j] = temp.est.att[j,3]
            YMAX[j] = temp.est.att[j,4]
            XVAR[j] = temp.time[j]
            NCOUNT[j] = temp.est.att[j,6]
          }
          data.toplot.main <- cbind.data.frame(YVAR,YMIN,YMAX,XVAR)
          if (length(ylim) != 0) {
            rect.length <- (ylim[2] - ylim[1]) / 5
            rect.min <- ylim[1]
          } else {
            rect.length <- (max(c(YMAX), na.rm = TRUE) - min(c(YMIN), na.rm = TRUE))/2
            rect.min <- min(c(YMIN), na.rm = TRUE) - rect.length
          }
        }

        temp.p <- ggplot(data.toplot.main)
        ## xlab and ylab
        temp.p <- temp.p + xlab(xlab) +  ylab(ylab)

        ## theme
        if (theme.bw == TRUE) {
          temp.p <- temp.p + theme_bw()
        }

        ## grid
        if (gridOff == TRUE) {
          temp.p <- temp.p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }

        # horizontal 0 line
        #temp.p <- temp.p + geom_hline(yintercetemp.pt = 0, colour = lcolor,size = lwidth)


        if(CI==FALSE){
          #temp.p <- temp.p + geom_point(aes(x=XVAR,y=YVAR,color='gray50',fill='gray50',alpha=1,size=1.2))
          temp.p <- temp.p + geom_point(aes(x=XVAR,y=YVAR),color='gray50',fill='gray50',alpha=1)
        } else {
          #temp.p <- temp.p + geom_pointrange(aes(x=XVAR,y=YVAR,ymin=YMIN,ymax=YMAX,color='gray50',fill='gray50',alpha=1,size=0.6))
          temp.p <- temp.p + geom_pointrange(aes(x=XVAR,y=YVAR,ymin=YMIN,ymax=YMAX),color='gray50',fill='gray50',alpha=1)
        }
        if(count==TRUE){

          T.start <- c()
          T.end <- c()
          ymin <- c()
          ymax <- c()
          T.gap <- (max(XVAR)-min(XVAR))/length(XVAR)
          for(j in c(1:length(XVAR))){
            T.start <- c(T.start,XVAR[j]-0.25*T.gap)
            T.end <- c(T.end,XVAR[j]+0.25*T.gap)
            ymin <- c(ymin, rect.min)
            ymax <- c(ymax, rect.min+rect.length*NCOUNT[j]/max(NCOUNT))
          }
          data.toplot <- cbind.data.frame(xmin=T.start,
                                          xmax=T.end,
                                          ymin=ymin,
                                          ymax=ymax)
          max.count.pos <- mean(XVAR[which.max(NCOUNT)])
          temp.p <- temp.p + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),data=data.toplot,fill='gray50',alpha=0.3,size=0.3,color='black')
          temp.p <- temp.p + annotate("text", x = max.count.pos - 0.02 * T.gap,
                            y = max(data.toplot$ymax) + 0.2 * rect.length,
                            label = max(NCOUNT), size = cex.text * 0.8, hjust = 0.5)
        }
        ## title
        if (is.null(main) == TRUE) {
          temp.title = paste0("Dynamic Effect when Moderator = ",x$Val.HTE[i])
          temp.p <- temp.p + ggtitle(temp.title) + theme(plot.title = element_text(hjust = 0.5))
        } else if (main!=""){
          temp.p <- temp.p + ggtitle(main)
        }
        p_collection[[i]] = temp.p

        if (is.null(save_path) == FALSE){
          filename <- file.path(save_path, paste0("plot_dynamic_",i,".png"))
          ggsave(filename,temp.p)
        }
    }
    p <- plot_grid(plotlist = p_collection, ncol = 1)
}
## ylim
# if (is.null(ylim) == FALSE) {
#   p <- p + coord_cartesian(ylim = ylim)
# }

# if(length(XVAR)<=10){
#   p <- p + scale_x_continuous(breaks=XVAR)
# } else {
#   p <- p + scale_x_continuous(labels=scaleFUN)
# }

# ## xlim
# if(is.null(xlim)){
#     if(is.na(d1[1,1])){
#         ## drop all periods before first non-missing
#         for(j in c(2:dim(d1)[1])){
#             if(!is.na(d1[j,1])){
#                 xlim <- c(XVAR[j],max(XVAR))
#                 break
#             }
#         }
#     }
# }
## xlim
# if (is.null(xlim) == FALSE) {
#   p <- p + coord_cartesian(xlim = xlim)
# }


# p <- p + geom_hline(yintercept = x$att.avg,color='red',size=0.8,linetype='dashed')
#
# p <- p + theme(legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = cex.legend),
#                legend.position = legend.pos,
#                legend.background = element_rect(fill="transparent",colour=NA),
#                axis.title=element_text(size=cex.lab),
#                axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
#                axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
#                axis.text = element_text(color="black", size=cex.axis),
#                #axis.text.x = element_text(size = cex.axis, angle = angle, hjust=x.h, vjust=x.v),
#                axis.text.x = element_text(size = cex.axis, angle = angle),
#                axis.text.y = element_text(size = cex.axis),
#                plot.title = element_text(size = cex.main, hjust = 0.5, face="bold", margin = margin(10, 0, 10, 0)))
return(p)
}
