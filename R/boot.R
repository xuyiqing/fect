###############################################
## Inference 
###############################################
fect.boot <- function(Y,
                      X,
                      D,
                      cl = NULL,
                      I,
                      II, 
                      T.on, 
                      T.off = NULL, 
                      T.on.carry = NULL,
                      T.on.balance = NULL,
                      balance.period = NULL,
                      method = "ife",
                      degree = 2,
                      sfe = NULL,
                      cfe = NULL,
                      ind.matrix = NULL,
                      knots = NULL,
                      criterion = "mspe",
                      CV,
                      k = 5,
                      cv.prop = 0.1, 
                      cv.treat = 0, 
                      cv.nobs = 1,
                      r = 0, 
                      r.end,
                      lambda = NULL,
                      nlambda = 10,
                      alpha = 0.05,
                      binary,
                      QR,
                      force,                      
                      hasRevs = 1,
                      tol,
                      norm.para,
                      placebo.period = NULL,
                      placeboTest = FALSE,
                      carryoverTest = FALSE,
                      carryover.period = NULL,
                      vartype = "bootstrap",
                      nboots,
                      parallel = TRUE,
                      cores = NULL,
                      group.level = NULL,
                      group = NULL,
                      dis = TRUE) {
    
    
    na.pos <- NULL
    TT <- dim(Y)[1]
    N <- dim(Y)[2]
    if (is.null(X) == FALSE) {
        p <- dim(X)[3]
    } else {
        p <- 0
    }


    if (hasRevs == 1) {
        ## D.fake : check reversals
        D.fake <- apply(D, 2, function(vec){cumsum(vec)})
        D.fake <- ifelse(D.fake > 0, 1, 0)
        D.fake[which(I == 0)] <- 0

        rev <- which(apply(D.fake == D, 2, sum) != TT)
        co <- which(apply(D, 2, sum) == 0)
        tr.all <- which(apply(D, 2, sum) > 0)
        tr <- tr.all[which(!tr.all %in% rev)]

        Nrev <- length(rev)
        Ntr <- length(tr)
        Nco <- length(co)
    } 
    else {
        ## treatement indicator
        tr <- which(apply(D, 2, sum) > 0)
        co <- which(apply(D, 2, sum) == 0)
        Ntr <- length(tr)
        Nco <- length(co)
        
    }

    
    ## estimation
    if (CV == 0) { 
        if(method == "gsynth"){
            out <- fect.gsynth(Y = Y, X = X, D = D, I = I, II = II, 
                           T.on = T.on, T.off = T.off, CV = 0, 
                           T.on.balance = T.on.balance,
                           balance.period = balance.period,
                           r = r, binary = binary, QR = QR,
                           force = force, hasRevs = hasRevs, 
                           tol = tol, boot = 0,
                           norm.para = norm.para, 
                           placebo.period = placebo.period,
                           placeboTest = placeboTest,
                           carryover.period = carryover.period,
                           carryoverTest = carryoverTest,
                           group.level = group.level, group = group)

        } else if (method == "ife") {
            out <- fect.fe(Y = Y, X = X, D = D, I = I, II = II, 
                           T.on = T.on, T.off = T.off, T.on.carry = T.on.carry,
                           T.on.balance = T.on.balance,
                           balance.period = balance.period,
                           r.cv = r, binary = binary, QR = QR,
                           force = force, hasRevs = hasRevs, 
                           tol = tol, boot = 0,
                           norm.para = norm.para, 
                           placebo.period = placebo.period,
                           placeboTest = placeboTest,
                           carryover.period = carryover.period,
                           carryoverTest = carryoverTest,
                           group.level = group.level, group = group)
        
        } else if (method == "mc") {
            out <- try(fect.mc(Y = Y, X = X, D = D, I = I, II = II,
                           T.on = T.on, T.off = T.off,  T.on.carry = T.on.carry,
                           T.on.balance = T.on.balance,
                           balance.period = balance.period,
                           lambda.cv = lambda, force = force, hasRevs = hasRevs, 
                           tol = tol, boot = 0,
                           norm.para = norm.para,
                           placebo.period = placebo.period,
                           placeboTest = placeboTest,
                           carryover.period = carryover.period,
                           carryoverTest = carryoverTest,
                           group.level = group.level, group = group), silent = TRUE)
            if ('try-error' %in% class(out)) {
                stop("\nCannot estimate using full data with MC algorithm.\n")
            }
        } else if (method %in% c("polynomial", "bspline","cfe")) {
            out <- try(fect.polynomial(Y = Y, D = D, X = X, I = I, 
                                   II = II, T.on = T.on,  T.on.carry = T.on.carry,
                                   T.on.balance = T.on.balance,
                                   balance.period = balance.period,
                                   T.off = T.off,
                                   method = method,degree = degree,
                                   knots = knots, force = force, 
                                   sfe = sfe, cfe = cfe,
                                   ind.matrix = ind.matrix,
                                   hasRevs = hasRevs,
                                   tol = tol, boot = 0, 
                                   placeboTest = placeboTest,
                                   placebo.period = placebo.period, 
                                   carryover.period = carryover.period,
                                   carryoverTest = carryoverTest,
                                   norm.para = norm.para,
                                   group.level = group.level, 
                                   group = group),silent = TRUE)
            #I.report <- out$I
            #II.report <- out$II
            if ('try-error' %in% class(out)) {
                stop("\nCannot estimate.\n")
            }
        }
    } 
    else {
        ## cross-valiadtion 
        if (binary == 0) {
            out <- fect.cv(Y = Y, X = X, D = D, I = I, II = II, 
                       T.on = T.on, T.off = T.off, T.on.carry = T.on.carry,
                       T.on.balance = T.on.balance,
                       balance.period = balance.period,
                       method = method, criterion = criterion,
                       k = k, r = r, r.end = r.end, 
                       nlambda = nlambda, lambda = lambda, 
                       force = force, hasRevs = hasRevs, 
                       tol = tol, norm.para = norm.para,
                       group.level = group.level, group = group,
                       cv.prop = cv.prop, cv.treat = cv.treat, 
                       cv.nobs = cv.nobs)

            method <- out$method
        } 
        else {
            out <- fect.binary.cv(Y = Y, X = X, D = D, 
                                  I = I, II = II, 
                                  T.on = T.on, T.off = T.off, 
                                  k = k, r = r, r.end = r.end, 
                                  QR = QR, force = force, 
                                  hasRevs = hasRevs, tol = tol,
                                  group.level = group.level, group = group)
            method <- "ife"
        }
        
    }
    
    
    
    ## output
    validX <- out$validX
    eff <- out$eff
    att.avg <- out$att.avg
    att.avg.unit <- out$att.avg.unit
    calendar.eff <- out$eff.calendar
    calendar.eff.fit <- out$eff.calendar.fit
    calendar.N <- out$N.calendar

    group.att <- out$group.att

    att <- out$att
    time.on <- out$time
    target.enp <- out$calendar.enp

    time.off <- NULL
    if (hasRevs == 1) {
        att.off <- out$att.off
        time.off <- out$time.off
    }
    carry.att <- carry.time <- NULL
    if (!is.null(T.on.carry)) {
      carry.att <- out$carry.att
      carry.time <- out$carry.time 
    }

    if(!is.null(balance.period)){
        balance.att <- out$balance.att
        balance.time <- out$balance.time
        balance.count <- out$balance.count
        balance.avg.att <- out$balance.avg.att
        if (!is.null(placebo.period) & placeboTest == TRUE) {
            balance.att.placebo <- out$balance.att.placebo
        }
    }

    eff.out <- out$eff
    fit.out <- out$Y.ct
    N_unit <- dim(out$res)[2]

    if(!is.null(group)){
        group.output.origin <- out$group.output
        group.output.name <- names(out$group.output)
        group.time.on <- list()
        group.time.off <- list()
        for(sub.name in group.output.name){
            group.time.on[[sub.name]] <- out$group.output[[sub.name]]$time.on
            if(hasRevs == 1){
                group.time.off[[sub.name]] <- out$group.output[[sub.name]]$time.off
            }
        }
    }

    if (p > 0) {
        beta <- out$beta
    } else {
        beta <-matrix(0,1,0)
    }

    if (is.null(cl)) {
        cl.unique <- NULL
    } else {
        cl.unique <- unique(cl)
    }

    if (vartype == "jackknife") {
        nboots <- N
    }
 
    ## bootstrapped estimates
    ## eff.boot <- array(0,dim = c(TT, Ntr, nboots))  ## to store results
    att.avg.boot <- matrix(0, 1, nboots)
    att.avg.unit.boot <- matrix(0, 1, nboots)
    att.boot <- matrix(0, length(time.on), nboots)
    att.count.boot <- matrix(0, length(time.on), nboots)
    beta.boot <- marginal.boot <- att.off.boot <- att.off.count.boot <- NULL
    calendar.eff.boot <- matrix(0,TT,nboots)
    calendar.eff.fit.boot <- matrix(0,TT,nboots)

    if (hasRevs == 1) {
        att.off.boot <- matrix(0, length(time.off), nboots) 
        att.off.count.boot <- matrix(0, length(time.off), nboots)   
    }
    if (!is.null(T.on.carry)) {
        carry.att.boot <- matrix(0, length(carry.att), nboots)
    }
    if (!is.null(balance.period)) {
        balance.att.boot <- matrix(0, length(balance.att), nboots) # dynamic att
        balance.count.boot <- matrix(0, length(balance.att), nboots) 
        balance.avg.att.boot <- matrix(0, 1, nboots)
        if (!is.null(placebo.period) & placeboTest == TRUE) {
            balance.att.placebo.boot <- matrix(0, 1, nboots)
        }
    }

    if (p > 0) {
        beta.boot <- matrix(0, p, nboots)
        if (binary == TRUE) {
            marginal.boot <- matrix(0, p, nboots)
        }
    }
    if (!is.null(placebo.period) & placeboTest == TRUE) {
        att.placebo.boot <- matrix(0, 1, nboots)
    }
    if (!is.null(carryover.period) & carryoverTest == TRUE) {
        att.carryover.boot <- matrix(0, 1, nboots)
    }

    group.att.boot <- NULL
    if (!is.null(group.att)) {
        group.att.boot <- matrix(0, length(group.att), nboots)
    }

    group.atts.boot <- NULL
    group.atts.off.boot <- NULL
    group.atts.placebo.boot <-NULL
    group.atts.carryover.boot <- NULL
    if(!is.null(group)){
        group.atts.boot <- list()
        group.atts.off.boot <- list()
        group.att.placebo.boot <- list()
        group.att.carryover.boot <- list()
        for(sub.name in group.output.name){
            subgroup.time.on <- group.time.on[[sub.name]]
            group.atts.boot[[sub.name]] <- matrix(0, length(subgroup.time.on), nboots)
            if(hasRevs == 1){
                subgroup.time.off <- group.time.off[[sub.name]]
                group.atts.off.boot[[sub.name]] <- matrix(0, length(subgroup.time.off), nboots)
            }
            if(placeboTest){
                group.att.placebo.boot[[sub.name]] <- matrix(0, 1, nboots)
            }
            if(carryoverTest){
                group.att.carryover.boot[[sub.name]] <- matrix(0, 1, nboots)
            }
        }   
    }

    if (dis) {
      if (vartype == "jackknife") {
          message("Jackknife estimates ... ")
      } else {
          message("Bootstrapping for uncertainties ... ")
      }
    }

    if (binary == TRUE & vartype == "parametric") {

        one.nonpara <- function(num = NULL) {
            Y.boot <- Y
            Y.fit <- out$Y.ct
            Y.fit[which(is.na(Y.fit))] <- 0
            Y.boot.new <- matrix(sapply(1:length(c(Y.fit)),function(i){rbinom(1,1,c(Y.fit)[i])}), TT, N)
            Y.boot[which(II == 1)] <- Y.boot.new[which(II == 1)]

            placebo.period.boot <- NULL
            carryover.period.boot <- NULL
            if (placeboTest == TRUE) {
                placebo.period.boot <- placebo.period
            }
            if (carryoverTest == TRUE) {
                carryover.period.boot <- carryover.period
            }

            boot <- try(fect.fe(Y = Y.boot, X = X, D = D,
                                    I = I, II = II, 
                                    T.on = T.on, T.off = T.off,  T.on.carry = T.on.carry,
                                    T.on.balance = T.on.balance,
                                    balance.period = balance.period,
                                    r.cv = out$r.cv, binary = binary,
                                    QR = QR, force = force,
                                    hasRevs = hasRevs, tol = tol, boot = 1,
                                    norm.para = norm.para,
                                    calendar.enp.seq = target.enp,
                                    time.on.seq = time.on, time.off.seq = time.off,
                                    time.on.carry.seq = carry.time,
                                    time.on.balance.seq = balance.time,
                                    placebo.period = placebo.period.boot, 
                                    placeboTest = placeboTest,
                                    carryoverTest = carryoverTest,
                                    carryover.period = carryover.period.boot,
                                    group.level = group.level,
                                    group = group), silent = TRUE)

            if ('try-error' %in% class(boot)) {
                boot0 <- list(att.avg = NA, att = NA, count = NA, 
                              beta = NA, att.off = NA, count.off = NA, eff.calendar = NA, 
                              eff.calendar.fit = NA,
                              att.placebo = NA, att.avg.unit = NA, att.carryover = NA,
                              group.att = NA, marginal = NA,carry.att = NA,balance.att = NA, 
                              balance.att.placebo = NA, balance.count = NA, 
                              balance.avg.att = NA, balance.time = NA,group.output = list())
                return(boot0)
            } else {
                return(boot)
            }

        }
    }
    else if(binary == FALSE & method %in% c("gsynth") & vartype == "parametric"){
        message("Parametric Bootstrap \n")
        sum.D <- colSums(out$D)
        id.tr <- which(sum.D>0)
        I.tr <- as.matrix(out$I[,id.tr])
        id.co <- which(sum.D==0)
        Nco <- length(id.co)
        Ntr <- length(id.tr)

        fit.out[which(out$I==0)] <- 0

        error.co <- out$res.full[,id.co]
        I.co <- out$I[,id.co]
        T0.ub <- apply(as.matrix(out$D[,id.tr] == 0), 2, sum) 
        T0.ub.min <- min(T0.ub)
        co.pre <- apply(as.matrix(I.co[1:T0.ub.min, ]), 2, sum)
        co.post <- apply(as.matrix(I.co[(max(T0.ub)+1):TT, ]), 2, sum)
        if (force %in% c(1, 3)) {
            valid.co <- id.co[(co.pre >= (out$r.cv + 1)) & (co.post >= 1)]
        } 
        else {
            valid.co <- id.co[(co.pre >= out$r.cv) & (co.post >= 1)]
        }

        draw.error <- function() {
            repeat {
                fake.tr <- sample(id.co, 1, replace = FALSE)
                if (fake.tr %in% valid.co) {
                    break
                }
            }
            
            id.co.rest <- id.co[which(!id.co %in% fake.tr)]
            repeat {
                id.co.pseudo <- sample(id.co.rest, Nco, replace = TRUE)
                if (sum(apply(as.matrix(out$I[, id.co.pseudo]), 1, sum) >= 1) == TT) {
                    break
                }
            }
                      
            id.pseudo <- c(rep(fake.tr, Ntr), id.co.pseudo)  ## Ntr + ...
            I.id.pseudo <- out$I[, id.pseudo] 
            II.id.pseudo <- out$II[,id.pseudo]
            ## obtain the prediction eror
            D.pseudo <- out$D[, c(id.tr, id.co.pseudo)]  ## fake.tr + control left
            Y.pseudo <- out$Y[, id.pseudo]
            T.on.pseudo <- T.on[,id.pseudo]
            
            X.pseudo <- NULL
            if (p > 0) {
                X.pseudo <- X[,id.pseudo,,drop = FALSE]
            }

            ## output
            #synth.out <- try(fect.gsynth(Y = Y.pseudo, X = X.pseudo, D = D.pseudo,
            #                             I = I.id.pseudo, II = II.id.pseudo,
            #                             force = force, r = out$r.cv, CV = 0,
            #                             tol = tol, norm.para = norm.para, boot = 1), silent = TRUE)

            synth.out <- try(fect.gsynth(Y = Y.pseudo, X = X.pseudo, D = D.pseudo,
                                     I = I.id.pseudo, II = II.id.pseudo,
                                     T.on = T.on.pseudo, hasRevs = hasRevs,
                                     force = force, r = out$r.cv, CV = 0,
                                     tol = tol, norm.para = norm.para, boot = 1), silent = TRUE)
   
            if ('try-error' %in% class(synth.out)) {
                return(matrix(NA, TT, Ntr))
            } 
            else {
                if ("eff" %in% names(synth.out)) {
                    if (is.null(norm.para)) {
                        output <- synth.out$eff.tr
                    } 
                    else {
                        output <- synth.out$eff.tr/norm.para[1]
                    }
                    
                    return(as.matrix(output)) ## TT * Ntr
                } 
                else {
                    return(matrix(NA, TT, Ntr))
                }
            }             
        }

        message("\rSimulating errors ...")
        if (parallel == TRUE) {
            error.tr <- foreach(j = 1:nboots,
                                .combine = function(...) abind(...,along=3),
                                .multicombine = TRUE,
                                .export = c("fect.gsynth","initialFit"),
                                .packages = c("fect","mvtnorm","fixest"),
                                .inorder = FALSE)  %dopar% {
                                    return(draw.error())
                                } 
        } 
        else {
            error.tr <- array(NA, dim = c(TT, Ntr, nboots))
            for (j in 1:nboots) {
                error.tr[,,j] <- draw.error()
                if (j %% 100 == 0) {
                    message(".")
                }
            }
        }



        if (0%in%I) {
            ## calculate vcov of ep_tr
            na.sum <- sapply(1:nboots, function(vec){sum(is.na(c(error.tr[,,vec])))})
            na.rm <- na.sum == TT * Ntr
            na.rm.count <- sum(na.rm)
            rm.pos <- which(na.rm == TRUE)

            if (na.rm.count > 0) {
                if (na.rm.count == nboots) {
                    stop("fail to simulate errors.\n")
                }
                error.tr <- error.tr[,,-rm.pos, drop = FALSE]
            }

            error.tr.adj <- array(NA, dim = c(TT, nboots - na.rm.count, Ntr))
            for(i in 1:Ntr){
                error.tr.adj[,,i] <- error.tr[,i,]
            }
            vcov_tr <- array(NA, dim = c(TT, TT, Ntr))
            for(i in 1:Ntr){
                vcov_tr[,,i] <- res.vcov(res = error.tr.adj[,,i],cov.ar = 0)
                vcov_tr[,,i][is.na(vcov_tr[,,i]) | is.nan(vcov_tr[,,i])] <- 0
            }  
            ## calculate vcov of e_co
            vcov_co <- res.vcov(res = error.co, cov.ar = 0)
            vcov_co[is.na(vcov_co) | is.nan(vcov_co)] <- 0
        }

        
        one.nonpara <- function(num = NULL){
            ## boostrap ID
            repeat {
                fake.co <- sample(id.co,Nco, replace=TRUE)
                if (sum(apply(as.matrix(I[,fake.co]), 1, sum) >= 1) == TT) {
                    break
                }
            }
            id.boot <- c(id.tr, fake.co)
                
            ## get the error for the treated and control
            error.tr.boot <- matrix(NA, TT, Ntr)
            if (0 %in% I) {        
                for (w in 1:Ntr) {
                    error.tr.boot[,w] <- t(rmvnorm(n = 1, rep(0, TT), vcov_tr[,,w], method = "svd"))
                } 
                error.tr.boot[which(I.tr == 0)] <- 0
                error.co.boot <- t(rmvnorm(n = Nco, rep(0, TT), vcov_co, method = "svd"))
                error.co.boot[which(as.matrix(I[,fake.co]) == 0)] <- 0    
            } 
            else {
                for (w in 1:Ntr) {
                    error.tr.boot[,w] <- error.tr[,w,sample(1:nboots,1,replace = TRUE)]
                }
                error.co.boot <- error.co[, sample(1:Nco, Nco, replace = TRUE)]   
            }

            Y.boot <- fit.out[,id.boot]
            Y.boot[,1:Ntr] <- as.matrix(Y.boot[,1:Ntr] + error.tr.boot)
            Y.boot[,(Ntr+1):length(id.boot)] <- Y.boot[,(Ntr+1):length(id.boot)] + error.co.boot 
            X.boot <- NULL
            if (p > 0) {
                X.boot <- X[,id.boot,,drop = FALSE] 
            }
            D.boot <- out$D[,id.boot] 
            I.boot <- out$I[,id.boot]
            II.boot <- out$II[,id.boot]
          
            synth.out <- try(fect.gsynth(Y = Y.boot, X = X.boot, D = D.boot,
                                         I = I.boot, II = II.boot,T.on = T.on[,id.boot], 
                                         T.on.balance = T.on.balance[,id.boot],
                                         balance.period = balance.period,
                                         hasRevs = hasRevs,
                                         force = force, r = out$r.cv, CV = 0, boot = 1,
                                         placeboTest = placeboTest,
                                         placebo.period = placebo.period, 
                                         carryover.period = carryover.period,
                                         carryoverTest = carryoverTest,
                                         calendar.enp.seq = target.enp,
                                         time.on.seq = time.on, time.off.seq = time.off,
                                         time.on.seq.group = group.time.on,
                                         time.off.seq.group = group.time.off,
                                         time.on.balance.seq = balance.time,
                                         norm.para = norm.para,tol = tol,
                                         group.level = group.level, group = group), silent = TRUE)

            if ('try-error' %in% class(synth.out)) {
                boot0 <- list(att.avg = NA, att = NA, count = NA, 
                                  beta = NA, att.off = NA, count.off = NA, eff.calendar = NA, 
                                  eff.calendar.fit = NA,
                                  att.placebo = NA, att.avg.unit = NA, att.carryover = NA,
                                  group.att = NA, marginal = NA,
                                  balance.att = NA, balance.att.placebo = NA, balance.count = NA, 
                                  balance.avg.att = NA, balance.time = NA,
                                  group.output = list())
                return(boot0)
            }
            else{
                return(synth.out)
            }
        }

    }
    else if(binary == FALSE & method %in% c("ife","mc","polynomial", "bspline","cfe") & vartype == 'parametric'){
        message("Parametric Bootstrap \n")
        sum.D <- colSums(out$D)
        tr <- which(sum.D>0)
        co <- which(sum.D==0)
        Nco <- length(co)
        Ntr <- length(tr)
        fit.out[which(out$I==0)] <- 0
        error.co <- out$res[,co]
        #error.tr <- out$eff[,tr]
        
        if (0%in%out$I) {
            vcov_co <- res.vcov(res = error.co, cov.ar = 0)
            vcov_co[is.na(vcov_co)|is.nan(vcov_co)] <- 0
            #vcov_tr <- res.vcov(res = error.tr, cov.ar = 0)
            #vcov_tr[is.na(vcov_tr)|is.nan(vcov_tr)] <- 0
        }
        
        one.nonpara <- function(num = NULL){
            error.id <- sample(1:Nco, N, replace = TRUE)    

            ## produce the new outcome data
            if (0%in%I) {
                error.boot <- t(rmvnorm(n=N,rep(0,TT),vcov_co,method="svd"))
                #error.boot.co <- t(rmvnorm(n=Nco,rep(0,TT),vcov_co,method="svd"))
                #error.boot.tr <- t(rmvnorm(n=Ntr,rep(0,TT),vcov_tr,method="svd"))
                Y.boot <- fit.out + out$eff + error.boot
                #Y.boot <- fit.out
                #Y.boot[,tr] <- Y.boot[,tr] +  error.boot.tr  
                #Y.boot[,co] <- Y.boot[,co] +  error.boot.co  
            } 
            else {
                Y.boot <- fit.out + out$eff + error.co[,error.id]
                #Y.boot <- fit.out
                #Y.boot[,tr] <- Y.boot[,tr] + error.tr[,error.id.tr]
                #Y.boot[,co] <- Y.boot[,co] + error.co[,error.id.co]
            }
            
            if (method == "ife") {
                boot <- try(fect.fe(Y = Y.boot, X = X, D = D, I = I, II = II, 
                            T.on = T.on, T.off = T.off,T.on.carry = T.on.carry,
                            T.on.balance = T.on.balance,
                            balance.period = balance.period,
                            r.cv = out$r.cv, binary = binary, QR = QR,
                            force = force, hasRevs = hasRevs, 
                            tol = tol, boot = 1,
                            norm.para = norm.para, 
                            placebo.period = placebo.period,
                            placeboTest = placeboTest,
                            carryover.period = carryover.period,
                            carryoverTest = carryoverTest,
                            group.level = group.level, group = group,
                            calendar.enp.seq = target.enp,
                            time.on.seq = time.on, time.off.seq = time.off,
                            time.on.carry.seq = carry.time,
                            time.on.balance.seq = balance.time,
                            time.on.seq.group = group.time.on,
                            time.off.seq.group = group.time.off),silent = TRUE)
            
            } 
            else if (method == "mc") {
                
                boot <- try(fect.mc(Y = Y.boot, X = X, D = D, I = I, II = II,
                            T.on = T.on, T.off = T.off, T.on.carry = T.on.carry,
                            T.on.balance = T.on.balance,
                            balance.period = balance.period,
                            lambda.cv = out$lambda.cv, force = force, hasRevs = hasRevs, 
                            tol = tol, boot = 1,
                            norm.para = norm.para,
                            placebo.period = placebo.period,
                            placeboTest = placeboTest,
                            carryover.period = carryover.period,
                            carryoverTest = carryoverTest,
                            group.level = group.level, group = group,
                            calendar.enp.seq = target.enp,
                            time.on.seq = time.on, time.off.seq = time.off,
                            time.on.carry.seq = carry.time,
                            time.on.balance.seq = balance.time,
                            time.on.seq.group = group.time.on,
                            time.off.seq.group = group.time.off),silent = TRUE)

            } 
            else if (method %in% c("polynomial", "bspline","cfe")) {
                boot <- try(fect.polynomial(Y = Y.boot, D = D, X = X, I = I, 
                                    II = II, T.on = T.on, 
                                    T.off = T.off,T.on.carry = T.on.carry,
                                    T.on.balance = T.on.balance,
                                    balance.period = balance.period,
                                    method = method,degree = degree,
                                    knots = knots, force = force,
                                    sfe = sfe, cfe = cfe,
                                    ind.matrix = ind.matrix, 
                                    hasRevs = hasRevs,
                                    tol = tol, boot = 1, 
                                    placeboTest = placeboTest,
                                    placebo.period = placebo.period, 
                                    carryover.period = carryover.period,
                                    carryoverTest = carryoverTest,
                                    norm.para = norm.para,
                                    group.level = group.level, group = group,
                                    calendar.enp.seq = target.enp,
                                    time.on.seq = time.on, time.off.seq = time.off,
                                    time.on.carry.seq = carry.time,
                                    time.on.balance.seq = balance.time,
                                    time.on.seq.group = group.time.on,
                                    time.off.seq.group = group.time.off), silent = TRUE)
            }

            if ('try-error' %in% class(boot)) {
                boot0 <- list(att.avg = NA, att = NA, count = NA, 
                                  beta = NA, att.off = NA, count.off = NA, eff.calendar = NA, 
                                  eff.calendar.fit = NA,
                                  att.placebo = NA, att.avg.unit = NA, att.carryover = NA,
                                  group.att = NA, marginal = NA,carry.att = NA, 
                                  balance.att = NA, balance.att.placebo = NA, balance.count = NA, 
                                  balance.avg.att = NA, balance.time = NA,
                                  group.output = list())
                return(boot0)
            }
            else{
                return(boot)
            }
        }
    } 
    else {
        
        one.nonpara <- function(num = NULL) { ## bootstrap
            if (is.null(num)) {
                if (is.null(cl)) {
                    if (hasRevs == 0) {
                        if (Nco > 0) {
                            repeat{
                                fake.co <- sample(co, Nco, replace=TRUE)
                                fake.tr <- sample(tr, Ntr, replace=TRUE)
                                boot.id <- c(fake.tr, fake.co)
                                if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
                                    break
                                }
                            }
                        } 
                        else {
                            repeat{
                                boot.id <- sample(tr, Ntr, replace=TRUE)
                                if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
                                    break
                                }
                            }
                        }
                    } 
                    else {
                        if (Ntr > 0) {
                            if (Nco > 0) {
                                repeat{
                                    fake.co <- sample(co, Nco, replace=TRUE)
                                    fake.tr <- sample(tr, Ntr, replace=TRUE)
                                    fake.rev <- sample(rev, Nrev, replace=TRUE)
                                    boot.id <- c(fake.rev, fake.tr, fake.co)
                                    if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
                                        break
                                    }
                                }
                            } 
                            else {
                                repeat{
                                    fake.tr <- sample(tr, Ntr, replace=TRUE)
                                    fake.rev <- sample(rev, Nrev, replace=TRUE)
                                    boot.id <- c(fake.rev, fake.tr)
                                    if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
                                        break
                                    }
                                }
                            }
                        } 
                        else {
                            if (Nco > 0) {
                                repeat{
                                    fake.co <- sample(co, Nco, replace=TRUE)
                                    fake.rev <- sample(rev, Nrev, replace=TRUE)
                                    boot.id <- c(fake.rev, fake.co)
                                    if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
                                        break
                                    }
                                }
                            } 
                            else {
                                repeat{
                                    boot.id <- sample(rev, Nrev, replace=TRUE)
                                    if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
                                        break
                                    }
                                }
                            }
                        }
                    }
                } 
                else {
                    cl.boot <- sample(cl.unique, length(cl.unique), replace = TRUE)
                    cl.boot.uni <- unique(cl.boot)
                    cl.boot.count <- as.numeric(table(cl.boot))
                    boot.id <- c()
                    for (kk in 1:length(cl.boot.uni)) {
                        boot.id <- c(boot.id, rep(which(cl == cl.boot.uni[kk]), cl.boot.count[kk]))
                    }
                }

                boot.group <- group[, boot.id]

            } 
            else { ## jackknife
                boot.group <- group[,-num]
                boot.id <- 1:N
                boot.id <- boot.id[-num]
            }
            
            X.boot <- X[,boot.id,,drop = FALSE]
            D.boot <- D[, boot.id]
            I.boot <- I[, boot.id]

            if(method=='cfe'){
                ind.matrix.boot <- list()
                for(ind.name in names(ind.matrix)){
                    ind.matrix.boot[[ind.name]] <- as.matrix(ind.matrix[[ind.name]][,boot.id])
                }
            }

            if (sum(c(D.boot) == 0) == 0 | sum(c(D.boot) == 1) == 0 | sum(c(I.boot) == 1) == 0) {
                boot0 <- list(att.avg = NA, att = NA, count = NA, 
                              beta = NA, att.off = NA, count.off = NA, eff.calendar = NA, 
                              eff.calendar.fit = NA,
                              att.placebo = NA, att.avg.unit = NA, att.carryover = NA,
                              group.att = list(),
                              balance.att = NA, balance.att.placebo = NA, balance.count = NA, 
                              balance.avg.att = NA, balance.time = NA)
                return(boot0)
            } 
            else {
                T.off.boot <- NULL
                if (hasRevs == TRUE) {
                    T.off.boot <- T.off[, boot.id]
                }
                placebo.period.boot <- NULL
                if (placeboTest == TRUE) {
                    placebo.period.boot <- placebo.period
                }
                carryover.period.boot <- NULL
                if(carryoverTest == TRUE){
                    carryover.period.boot <- carryover.period
                }

                if(method == "gsynth") {
                    boot <- try(fect.gsynth(Y = Y[, boot.id], X = X.boot, D = D.boot,
                                    I = I.boot, II = II[, boot.id], 
                                    T.on = T.on[, boot.id], T.off = T.off.boot, CV = 0,
                                    T.on.balance = T.on.balance[, boot.id],
                                    balance.period = balance.period,
                                    r = out$r.cv, binary = binary,
                                    QR = QR, force = force,
                                    hasRevs = hasRevs, tol = tol, boot = 1,
                                    norm.para = norm.para,
                                    calendar.enp.seq = target.enp,
                                    time.on.seq = time.on, time.off.seq = time.off,
                                    placebo.period = placebo.period.boot, 
                                    placeboTest = placeboTest,
                                    time.on.balance.seq = balance.time,
                                    carryoverTest = carryoverTest,
                                    carryover.period = carryover.period.boot,
                                    group.level = group.level,
                                    group = boot.group,
                                    time.on.seq.group = group.time.on,
                                    time.off.seq.group = group.time.off) )           
                }
                else if (method == "ife") {
                    boot <- try(fect.fe(Y = Y[, boot.id], X = X.boot, D = D.boot,
                                    I = I.boot, II = II[, boot.id], 
                                    T.on = T.on[, boot.id], T.off = T.off.boot, 
                                    T.on.carry = T.on.carry[, boot.id],
                                    T.on.balance = T.on.balance[, boot.id],
                                    balance.period = balance.period,
                                    r.cv = out$r.cv, binary = binary,
                                    QR = QR, force = force,
                                    hasRevs = hasRevs, tol = tol, boot = 1,
                                    norm.para = norm.para,
                                    calendar.enp.seq = target.enp,
                                    time.on.seq = time.on, time.off.seq = time.off,
                                    time.on.carry.seq = carry.time,
                                    time.on.balance.seq = balance.time,
                                    placebo.period = placebo.period.boot, 
                                    placeboTest = placeboTest,
                                    carryoverTest = carryoverTest,
                                    carryover.period = carryover.period.boot,
                                    group.level = group.level,
                                    group = boot.group,
                                    time.on.seq.group = group.time.on,
                                    time.off.seq.group = group.time.off), silent = TRUE)
                } else if (method == "mc") {
                    boot <- try(fect.mc(Y = Y[,boot.id], X = X.boot, D = D[,boot.id],
                                    I = I[,boot.id], II = II[,boot.id],
                                    T.on = T.on[,boot.id], T.off = T.off.boot,  
                                    T.on.carry = T.on.carry[, boot.id],
                                    T.on.balance = T.on.balance[, boot.id],
                                    balance.period = balance.period,
                                    lambda.cv = out$lambda.cv, force = force, 
                                    hasF = out$validF, hasRevs = hasRevs, 
                                    tol = tol, boot = 1,
                                    norm.para = norm.para,calendar.enp.seq = target.enp,
                                    time.on.seq = time.on, time.off.seq = time.off,
                                    time.on.carry.seq = carry.time,
                                    time.on.balance.seq = balance.time,
                                    placebo.period = placebo.period.boot, 
                                    placeboTest = placeboTest,
                                    carryoverTest = carryoverTest,
                                    carryover.period = carryover.period.boot,
                                    group.level = group.level,
                                    group = boot.group,
                                    time.on.seq.group = group.time.on,
                                    time.off.seq.group = group.time.off), silent = TRUE)

                } else if (method %in% c("polynomial", "bspline", "cfe")) {
                    
                    boot <- try(fect.polynomial(Y = Y[,boot.id], X = X.boot, 
                                                    D = D[,boot.id],
                                                    I = I[,boot.id], II = II[,boot.id],
                                                    T.on = T.on[,boot.id], T.off = T.off.boot,  
                                                    T.on.carry = T.on.carry[, boot.id],
                                                    T.on.balance = T.on.balance[, boot.id],
                                                    balance.period = balance.period,
                                                    method = method, degree = degree, 
                                                    sfe = sfe, cfe = cfe,
                                                    ind.matrix = ind.matrix.boot,
                                                    knots = knots,
                                                    force = force, hasRevs = hasRevs,
                                                    tol = tol,boot = 1,
                                                    norm.para = norm.para, time.on.seq = time.on, 
                                                    calendar.enp.seq = target.enp,
                                                    time.off.seq = time.off,
                                                    time.on.carry.seq = carry.time,
                                                    time.on.balance.seq = balance.time,
                                                    placebo.period = placebo.period.boot, 
                                                    carryoverTest = carryoverTest,
                                                    carryover.period = carryover.period.boot,
                                                    placeboTest = placeboTest,
                                                    group.level = group.level,
                                                    group = boot.group,
                                                    time.on.seq.group = group.time.on,
                                                    time.off.seq.group = group.time.off),silent = TRUE)
                                            
                                 
                }

                if ('try-error' %in% class(boot)) {
                    boot0 <- list(att.avg = NA, att = NA, count = NA, 
                                  beta = NA, att.off = NA, count.off = NA, eff.calendar = NA, 
                                  eff.calendar.fit = NA,
                                  att.placebo = NA, att.avg.unit = NA, att.carryover = NA,
                                  group.att = NA, marginal = NA,carry.att = NA,
                                  group.output = list(),
                                  balance.att = NA, balance.att.placebo = NA, balance.count = NA, 
                                  balance.avg.att = NA, balance.time = NA)
                    return(boot0)
                } else {
                    return(boot)
                }
            }            
        }
    } 
            
    
    ## jack.seq <- sample(1:N, N, replace = FALSE)
    boot.seq <- NULL
    if (vartype == "jackknife") {
        ## nboots <- min(N, nboots)
        ## boot.seq <- jack.seq[1:nboots]
        boot.seq <- 1:N 
    }


    ## computing
    if (parallel == TRUE) { 
        boot.out <- foreach(j=1:nboots, 
                            .inorder = FALSE,
                            .export = c("fect.fe", "fect.mc", "fect.polynomial", "get_term","fect.gsynth","initialFit"),
                            .packages = c("fect","mvtnorm","fixest")
                            ) %dopar% {
                                return(one.nonpara(boot.seq[j]))
                            }

        for (j in 1:nboots) { 
            att.avg.boot[,j] <- boot.out[[j]]$att.avg
            att.avg.unit.boot[, j] <- boot.out[[j]]$att.avg.unit
            att.boot[,j] <- boot.out[[j]]$att
            att.count.boot[,j] <- boot.out[[j]]$count
            
            calendar.eff.boot[,j] <- boot.out[[j]]$eff.calendar
            calendar.eff.fit.boot[,j] <- boot.out[[j]]$eff.calendar.fit
            if (p > 0) {
                beta.boot[,j] <- boot.out[[j]]$beta
                if (binary == TRUE) {
                    marginal.boot[, j] <- boot.out[[j]]$marginal
                }
            }
            if (hasRevs == 1) {
                att.off.boot[,j] <- boot.out[[j]]$att.off
                att.off.count.boot[,j] <- boot.out[[j]]$count.off 
            }
            if (!is.null(T.on.carry)) {
                carry.att.boot[,j] <- boot.out[[j]]$carry.att
            }
            if(!is.null(balance.period)){
                balance.att.boot[,j] <- boot.out[[j]]$balance.att
                balance.count.boot[,j] <- boot.out[[j]]$balance.count
                balance.avg.att.boot[,j] <- boot.out[[j]]$balance.avg.att
                if (!is.null(placebo.period) & placeboTest == TRUE) {
                    balance.att.placebo.boot[,j] <- boot.out[[j]]$balance.att.placebo
                }
            }

            if (!is.null(placebo.period) & placeboTest == TRUE) {
                att.placebo.boot[,j] <- boot.out[[j]]$att.placebo
            }
            if (!is.null(carryover.period) & carryoverTest == TRUE) {
                att.carryover.boot[,j] <- boot.out[[j]]$att.carryover
            }
            if (!is.null(group)) {
                group.att.boot[,j] <- boot.out[[j]]$group.att  
                for(sub.name in group.output.name){
                    if(is.null(boot.out[[j]]$group.output[[sub.name]]$att.on)){
                        group.atts.boot[[sub.name]][,j] <- NA
                    }else{
                        group.atts.boot[[sub.name]][,j] <- boot.out[[j]]$group.output[[sub.name]]$att.on
                    }
                    if(hasRevs == 1){
                        if(is.null(boot.out[[j]]$group.output[[sub.name]]$att.off)){
                            group.atts.off.boot[[sub.name]][,j] <- NA
                        }else{
                           group.atts.off.boot[[sub.name]][,j] <- boot.out[[j]]$group.output[[sub.name]]$att.off
                        }
                    }
                    if(placeboTest){
                        if(is.null(boot.out[[j]]$group.output[[sub.name]]$att.placebo)){
                            group.att.placebo.boot[[sub.name]][,j] <- NA
                        }else{
                           group.att.placebo.boot[[sub.name]][,j] <- boot.out[[j]]$group.output[[sub.name]]$att.placebo
                        }
                    }
                    if(carryoverTest){
                        if(is.null(boot.out[[j]]$group.output[[sub.name]]$att.carryover)){
                            group.att.carryover.boot[[sub.name]][,j] <- NA
                        }else{
                           group.att.carryover.boot[[sub.name]][,j] <- boot.out[[j]]$group.output[[sub.name]]$att.carryover
                        }
                    }
                }   
            }
        } 
    } 
    else {
        for (j in 1:nboots) { 
            boot <- one.nonpara(boot.seq[j]) 
            att.avg.boot[,j] <- boot$att.avg
            att.avg.unit.boot[,j] <- boot$att.avg.unit
            att.boot[,j] <- boot$att
            att.count.boot[,j] <- boot$count

            calendar.eff.boot[,j] <- boot$eff.calendar
            calendar.eff.fit.boot[,j] <- boot$eff.calendar.fit
            if (p > 0) {
                beta.boot[,j] <- boot$beta
                if (binary == TRUE) {
                    marginal.boot[,j] <- boot$marginal      
                }
            }
            if (hasRevs == 1) {
                att.off.boot[,j] <- boot$att.off
                att.off.count.boot[,j] <- boot$count.off 
            }
            if (!is.null(T.on.carry)) {
                carry.att.boot[,j] <- boot$carry.att
            }
            if(!is.null(balance.period)){
                balance.att.boot[,j] <- boot$balance.att
                balance.count.boot[,j] <- boot$balance.count
                balance.avg.att.boot[,j] <- boot$balance.avg.att
                if (!is.null(placebo.period) & placeboTest == TRUE) {
                    balance.att.placebo.boot[,j] <- boot$balance.att.placebo
                }
            }
            if (!is.null(placebo.period) & placeboTest == TRUE) {
                att.placebo.boot[,j] <- boot$att.placebo
            }
            if (!is.null(carryover.period) & carryoverTest == TRUE) {
                att.carryover.boot[,j] <- boot$att.carryover
            }
            if (!is.null(group)) {
                group.att.boot[,j] <- boot$group.att  
                for(sub.name in group.output.name){
                    if(is.null(boot$group.output[[sub.name]]$att.on)){
                        group.atts.boot[[sub.name]][,j] <- NA
                    }else{
                        group.atts.boot[[sub.name]][,j] <- boot$group.output[[sub.name]]$att.on
                    }
                    if(hasRevs == 1){
                        if(is.null(boot$group.output[[sub.name]]$att.off)){
                            group.atts.off.boot[[sub.name]][,j] <- NA
                        }else{
                           group.atts.off.boot[[sub.name]][,j] <- boot$group.output[[sub.name]]$att.off
                        }
                    }
                    if(placeboTest){
                        if(is.null(boot$group.output[[sub.name]]$att.placebo)){
                            group.att.placebo.boot[[sub.name]][,j] <- NA
                        }else{
                           group.att.placebo.boot[[sub.name]][,j] <- boot$group.output[[sub.name]]$att.placebo
                        }
                    }
                    if(carryoverTest){
                        if(is.null(boot$group.output[[sub.name]]$att.carryover)){
                            group.att.carryover.boot[[sub.name]][,j] <- NA
                        }else{
                           group.att.carryover.boot[[sub.name]][,j] <- boot$group.output[[sub.name]]$att.carryover
                        }
                    }
                }   
            }
            ## report progress
            if (j%%100 == 0)  {
                message(".")   
            }  
        }  
    } 
    ## end of bootstrapping

    ## remove failure bootstrap
    ## alternative condition? max(apply(is.na(att.boot),2,sum)) == dim(att.boot)[1]
    if (sum(is.na(c(att.avg.boot))) > 0) {
        boot.rm <- which(is.na(c(att.avg.boot)))
        att.avg.boot <- t(as.matrix(att.avg.boot[,-boot.rm]))
        att.avg.unit.boot <- t(as.matrix(att.avg.unit.boot[,-boot.rm]))
        att.boot <- as.matrix(att.boot[,-boot.rm])
        att.count.boot <- as.matrix(att.count.boot[,-boot.rm])
        calendar.eff.boot <- as.matrix(calendar.eff.boot[,-boot.rm])
        calendar.eff.fit.boot <- as.matrix(calendar.eff.fit.boot[,-boot.rm])
        if (p > 0) {
            beta.boot <- as.matrix(beta.boot[,-boot.rm])
            if (dim(beta.boot)[2] == 1) {
                beta.boot <- t(beta.boot)
            }
        }
        if (hasRevs == 1) {
            att.off.boot <- as.matrix(att.off.boot[,-boot.rm])
            att.off.count.boot <- as.matrix(att.off.count.boot[,-boot.rm])
        }
        if (!is.null(T.on.carry)) {
            carry.att.boot <- as.matrix(carry.att.boot[, -boot.rm])
        }
        if(!is.null(balance.period)){
            balance.att.boot <- as.matrix(balance.att.boot[, -boot.rm])
            balance.count.boot <- as.matrix(balance.count.boot[, -boot.rm])
            balance.avg.att.boot <- t(as.matrix(balance.avg.att.boot[,-boot.rm]))
            if (!is.null(placebo.period) & placeboTest == TRUE) {
                balance.att.placebo.boot <- t(as.matrix(balance.att.placebo.boot[,-boot.rm]))
            }
        }
        if (!is.null(placebo.period) & placeboTest == TRUE) {
            att.placebo.boot <- t(as.matrix(att.placebo.boot[,-boot.rm]))
        }
        if (!is.null(carryover.period) & carryoverTest == TRUE) {
            att.carryover.boot <- t(as.matrix(att.carryover.boot[,-boot.rm]))
        }
        if (!is.null(group)) {
            if (dim(group.att.boot)[1] == 1) {
                group.att.boot <- t(as.matrix(group.att.boot[, -boot.rm]))
            } else {
                group.att.boot <- as.matrix(group.att.boot[, -boot.rm])
            }

            for(sub.name in group.output.name){
                group.atts.boot[[sub.name]] <- as.matrix(group.atts.boot[[sub.name]][,-boot.rm])
                if(hasRevs == 1){
                    group.atts.off.boot[[sub.name]] <- as.matrix(group.atts.off.boot[[sub.name]][,-boot.rm])
                }
                if(placeboTest){
                    group.att.placebo.boot[[sub.name]] <- t(as.matrix(group.att.placebo.boot[[sub.name]][,-boot.rm]))
                }
                if(carryoverTest){
                    group.att.carryover.boot[[sub.name]] <- t(as.matrix(group.att.carryover.boot[[sub.name]][,-boot.rm]))
                }
            }  
        }
    }
    if (dis) {
        message(dim(att.boot)[2], " runs\n", sep = "")
    }


    ####################################
    ## Variance and CIs
    ####################################

    ## function to get two-sided p-values
    get.pvalue <- function(vec) {
        if (NaN%in%vec|NA%in%vec) {
            nan.pos <- is.nan(vec)
            na.pos <- is.na(vec)
            pos <- c(which(nan.pos),which(na.pos))
            vec.a <- vec[-pos]
            a <- sum(vec.a >= 0)/(length(vec)-sum(nan.pos|na.pos)) * 2
            b <- sum(vec.a <= 0)/(length(vec)-sum(nan.pos|na.pos)) * 2  
        } else {
            a <- sum(vec >= 0)/length(vec) * 2
            b <- sum(vec <= 0)/length(vec) * 2  
        }
        return(min(as.numeric(min(a, b)),1))
    }

    ## ATT estimates
    if (vartype == "jackknife") {
        
        att.j <- jackknifed(att, att.boot, alpha)
        est.att <- cbind(att, att.j$se, att.j$CI.l, att.j$CI.u, att.j$P, out$count)
        colnames(est.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                                  "p.value", "count")
        rownames(est.att) <- out$time

        att.bound <- cbind(att + qnorm(alpha)*att.j$se, att + qnorm(1 - alpha)*att.j$se)
        colnames(att.bound) <- c("CI.lower", "CI.upper")
        rownames(att.bound) <- out$time

        eff.calendar.j <- jackknifed(calendar.eff, calendar.eff.boot, alpha)
        est.eff.calendar <- cbind(calendar.eff, eff.calendar.j$se, eff.calendar.j$CI.l, eff.calendar.j$CI.u, eff.calendar.j$P, calendar.N)
        colnames(est.eff.calendar) <- c("ATT-calendar", "S.E.", "CI.lower", "CI.upper","p.value", "count")

        eff.calendar.fit.j <- jackknifed(calendar.eff.fit, calendar.eff.fit.boot, alpha)
        est.eff.calendar.fit <- cbind(calendar.eff.fit, eff.calendar.fit.j$se, eff.calendar.fit.j$CI.l, eff.calendar.fit.j$CI.u, eff.calendar.fit.j$P, calendar.N)
        colnames(est.eff.calendar.fit) <- c("ATT-calendar Fitted", "S.E.", "CI.lower", "CI.upper","p.value", "count")

        if (hasRevs == 1) {
            att.off.j <- jackknifed(att.off, att.off.boot, alpha)
            est.att.off <- cbind(att.off, att.off.j$se, att.off.j$CI.l, att.off.j$CI.u, att.off.j$P, out$count.off)
            colnames(est.att.off) <- c("ATT.OFF", "S.E.", "CI.lower", "CI.upper",
                                      "p.value", "count")
            rownames(est.att.off) <- out$time.off

            att.off.bound <- cbind(att.off + qnorm(alpha)*att.off.j$se, att.off + qnorm(1 - alpha)*att.off.j$se)
            colnames(att.off.bound) <- c("CI.lower", "CI.upper")
            rownames(att.off.bound) <- out$time.off
        }

        if (!is.null(T.on.carry)) {
            carry.att.j <- jackknifed(carry.att, carry.att.boot, alpha)
            est.carry.att <- cbind(carry.att, carry.att.j$se, 
                             carry.att.j$CI.l, carry.att.j$CI.u, carry.att.j$P)

            colnames(est.carry.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                                      "p.value")
            rownames(est.carry.att) <- carry.time
        }

        if(!is.null(balance.period)){
            balance.att.j <- jackknifed(balance.att, balance.att.boot, alpha)
            est.balance.att <- cbind(balance.att, balance.att.j$se, balance.att.j$CI.l, balance.att.j$CI.u, balance.att.j$P, out$balance.count)
            colnames(est.balance.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                                           "p.value", "count")
            rownames(est.balance.att) <- out$balance.time
            
            balance.avg.att.j <- jackknifed(balance.avg.att, balance.avg.att.boot, alpha)
            est.balance.avg <- t(as.matrix(c(balance.avg.att, balance.avg.att.j$se, balance.avg.att.j$CI.l, balance.avg.att.j$CI.u, balance.avg.att.j$P)))
            colnames(est.balance.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")

            balance.att.bound <- cbind(balance.att + qnorm(alpha)*balance.att.j$se, balance.att + qnorm(1 - alpha)*balance.att.j$se)
            colnames(balance.att.bound) <- c("CI.lower", "CI.upper")
            rownames(balance.att.bound) <- out$balance.time

            if (!is.null(placebo.period) & placeboTest == TRUE) {
                balance.att.placebo.j <- jackknifed(balance.att.placebo, balance.att.placebo.boot, alpha)
                est.balance.placebo <- t(as.matrix(c(balance.att.placebo, balance.att.placebo.j$se, balance.att.placebo.j$CI.l, balance.att.placebo.j$CI.u, balance.att.placebo.j$P)))
                colnames(est.balance.placebo) <- c("ATT.placebo", "S.E.", "CI.lower", "CI.upper", "p.value")
            } 

        }

        ## average (over time) ATT
        att.avg.j <- jackknifed(att.avg, att.avg.boot, alpha)
        est.avg <- t(as.matrix(c(att.avg, att.avg.j$se, att.avg.j$CI.l, att.avg.j$CI.u, att.avg.j$P)))
        colnames(est.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")

        att.avg.unit.j <- jackknifed(att.avg.unit, att.avg.unit.boot, alpha)
        est.avg.unit <- t(as.matrix(c(att.avg.unit, att.avg.unit.j$se, att.avg.unit.j$CI.l, att.avg.unit.j$CI.u, att.avg.unit.j$P)))
        colnames(est.avg.unit) <- c("ATT.avg.unit", "S.E.", "CI.lower", "CI.upper", "p.value")

        ## regression coefficents
        if (p > 0) {
            beta.j <- jackknifed(beta, beta.boot, alpha)
            est.beta <- cbind(beta, beta.j$se, beta.j$CI.l, beta.j$CI.u, beta.j$P)
            colnames(est.beta)<-c("beta", "S.E.", "CI.lower", "CI.upper", "p.value")

            if (binary == TRUE) {
                marginal.j <- jackknifed(out$marginal, marginal.boot, alpha)
                est.marginal <- cbind(out$marginal, marginal.j$se, marginal.j$CI.l, marginal.j$CI.u, marginal.j$P)
                colnames(est.marginal)<-c("marginal", "S.E.", "CI.lower", "CI.upper", "p.value")
            }
        }

        ## placebo test
        if (!is.null(placebo.period) & placeboTest == TRUE) {
            att.placebo <- out$att.placebo
            att.placebo.j <- jackknifed(att.placebo, att.placebo.boot, alpha)
            att.placebo.bound <- c(att.placebo + qnorm(alpha)*att.placebo.j$se, 
                                   att.placebo + qnorm(1 - alpha)*att.placebo.j$se)
            est.placebo <- t(as.matrix(c(att.placebo, att.placebo.j$se, 
                                         att.placebo.j$CI.l, att.placebo.j$CI.u, 
                                         att.placebo.j$P,
                                         att.placebo.bound)))
            colnames(est.placebo) <- c("ATT.placebo", "S.E.", 
                                       "CI.lower", "CI.upper", 
                                       "p.value", "CI.lower(90%)","CI.upper(90%)")
        }

        ## carryover test
        if (!is.null(carryover.period) & carryoverTest == TRUE) {
            att.carryover <- out$att.carryover
            att.carryover.j <- jackknifed(att.carryover, att.carryover.boot, alpha)
            att.carryover.bound <- c(att.carryover + qnorm(alpha)*att.carryover.j$se, 
                                   att.carryover + qnorm(1 - alpha)*att.carryover.j$se)
            
            est.carryover <- t(as.matrix(c(att.carryover, att.carryover.j$se, 
                                           att.carryover.j$CI.l, att.carryover.j$CI.u, 
                                           att.carryover.j$P,
                                           att.carryover.bound)))
            colnames(est.carryover) <- c("ATT.carryover", "S.E.", 
                                         "CI.lower", "CI.upper", 
                                         "p.value", "CI.lower(90%)","CI.upper(90%)")
        }

        ## cohort effect
        est.group.out <- NULL
        if (!is.null(group)) {
            group.att.j <- jackknifed(group.att, group.att.boot, alpha)
            est.group.att <- cbind(group.att, group.att.j$se, group.att.j$CI.l, group.att.j$CI.u, group.att.j$P)
            colnames(est.group.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                                         "p.value")
            
            est.group.out <- list()
            for(sub.name in group.output.name){
                subgroup.atts <- group.output.origin[[sub.name]]$att.on
                subgroup.atts.boot <- group.atts.boot[[sub.name]]
                subgroup.est.att <- NULL
                subgroup.att.bound <- NULL
                
                if(dim(subgroup.atts.boot)[1]>0){
                    subgroup.att.j <- jackknifed(subgroup.atts, subgroup.atts.boot, alpha)
                    subgroup.est.att <- cbind(subgroup.atts, subgroup.att.j$se, subgroup.att.j$CI.l, 
                                            subgroup.att.j$CI.u, subgroup.att.j$P, 
                                            group.output.origin[[sub.name]]$count.on)
                    colnames(subgroup.est.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                                        "p.value", "count")
                    rownames(subgroup.est.att) <- group.output.origin[[sub.name]]$time.on
                
                    subgroup.att.bound <- cbind(subgroup.atts + qnorm(alpha)*subgroup.att.j$se, 
                                                subgroup.atts + qnorm(1 - alpha)*subgroup.att.j$se)
                    colnames(subgroup.att.bound) <- c("CI.lower", "CI.upper")
                    rownames(subgroup.att.bound) <- group.output.origin[[sub.name]]$time.on
                }
                
                subgroup.est.att.off <- NULL
                subgroup.att.off.bound <- NULL
                if(hasRevs == 1){
                    subgroup.atts.off <- group.output.origin[[sub.name]]$att.off
                    subgroup.atts.off.boot <- group.atts.off.boot[[sub.name]]
                    if(dim(subgroup.atts.off.boot)[1]>0){
                        subgroup.att.off.j <- jackknifed(subgroup.atts.off, subgroup.atts.off.boot, alpha)
                        subgroup.est.att.off <- cbind(subgroup.atts.off, subgroup.att.off.j$se, subgroup.att.off.j$CI.l, 
                                                subgroup.att.off.j$CI.u, subgroup.att.off.j$P, 
                                                group.output.origin[[sub.name]]$count.off)
                        colnames(subgroup.est.att.off) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                                                            "p.value", "count")
                        rownames(subgroup.est.att.off) <- group.output.origin[[sub.name]]$time.off
                    
                        subgroup.att.off.bound <- cbind(subgroup.atts.off + qnorm(alpha)*subgroup.att.off.j$se, 
                                                    subgroup.atts.off + qnorm(1 - alpha)*subgroup.att.off.j$se)
                        colnames(subgroup.att.off.bound) <- c("CI.lower", "CI.upper")
                        rownames(subgroup.att.off.bound) <- group.output.origin[[sub.name]]$time.off                      
                    }
                }

                subgroup.est.placebo <- NULL
                if(placeboTest){
                    subgroup.att.placebo <- group.output.origin[[sub.name]]$att.placebo
                    if(length(subgroup.att.placebo)>0){
                        subgroup.att.placebo.j <- jackknifed(subgroup.att.placebo, group.att.placebo.boot[[sub.name]], alpha)
                        att.placebo.bound <- c(subgroup.att.placebo + qnorm(alpha)*subgroup.att.placebo.j$se, 
                                               subgroup.att.placebo + qnorm(1 - alpha)*subgroup.att.placebo.j$se)
            
                        subgroup.est.placebo <- t(as.matrix(c(subgroup.att.placebo, 
                                                            subgroup.att.placebo.j$se, 
                                                            subgroup.att.placebo.j$CI.l, 
                                                            subgroup.att.placebo.j$CI.u, 
                                                            subgroup.att.placebo.j$P,
                                                            att.placebo.bound)))
                        colnames(subgroup.est.placebo) <- c("ATT.placebo", "S.E.", 
                                                            "CI.lower", "CI.upper", "p.value",
                                                            "CI.lower(90%)","CI.upper(90%)")
                                            
                    }
                }

                subgroup.est.carryover <- NULL
                if(carryoverTest){
                    subgroup.att.carryover <- group.output.origin[[sub.name]]$att.carryover
                    if(length(subgroup.att.carryover)>0){
                        subgroup.att.carryover.j <- jackknifed(subgroup.att.carryover, group.att.carryover.boot[[sub.name]], alpha)
                        att.carryover.bound <- c(subgroup.att.carryover + qnorm(alpha)*subgroup.att.carryover.j$se, 
                                                 subgroup.att.carryover + qnorm(1 - alpha)*subgroup.att.carryover.j$se)
                                    
                        subgroup.est.carryover <- t(as.matrix(c(subgroup.att.carryover, 
                                                            subgroup.att.carryover.j$se, 
                                                            subgroup.att.carryover.j$CI.l, 
                                                            subgroup.att.carryover.j$CI.u, 
                                                            subgroup.att.carryover.j$P,
                                                            att.carryover.bound)))
                        colnames(subgroup.est.carryover) <- c("ATT.carryover", "S.E.", 
                                                              "CI.lower", "CI.upper", "p.value",
                                                              "CI.lower(90%)","CI.upper(90%)")
                                            
                    }
                }

                est.group.out[[sub.name]] <- list(att.on = subgroup.est.att,
                                                  att.on.bound = subgroup.att.bound,
                                                  att.on.boot = group.atts.boot[[sub.name]],
                                                  att.off = subgroup.est.att.off,
                                                  att.off.bound = subgroup.att.off.bound,
                                                  att.off.boot = group.atts.off.boot[[sub.name]],
                                                  att.placebo = subgroup.est.placebo,
                                                  att.carryover = subgroup.est.carryover)
            }
        }
    } 
    else {

        se.att <- apply(att.boot, 1, function(vec) sd(vec, na.rm=TRUE))
        CI.att <- cbind(att - se.att * qnorm(1-alpha/2), att + se.att * qnorm(1-alpha/2)) # normal approximation
        pvalue.att <- (1-pnorm(abs(att/se.att)))*2
        est.att <- cbind(att, se.att, CI.att, pvalue.att, out$count)
        colnames(est.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                                  "p.value", "count")
        rownames(est.att) <- out$time
        
        # for equivalence test
        att.bound <- cbind(att - se.att * qnorm(1-alpha), att + se.att * qnorm(1-alpha)) # one-sided
        colnames(att.bound) <- c("CI.lower", "CI.upper")
        rownames(att.bound) <- out$time


        
        if (hasRevs == 1) {
            CI.att.off <- t(apply(att.off.boot, 1, function(vec) quantile(vec,c(alpha/2, 1 - alpha/2), na.rm=TRUE)))
            se.att.off <- apply(att.off.boot, 1, function(vec) sd(vec, na.rm=TRUE))
            pvalue.att.off <- apply(att.off.boot, 1, get.pvalue)

            est.att.off <- cbind(att.off, se.att.off, CI.att.off, pvalue.att.off, out$count.off)
            colnames(est.att.off) <- c("ATT.OFF", "S.E.", "CI.lower", "CI.upper",
                                       "p.value", "count.off")
            rownames(est.att.off) <- out$time.off
            #T0.off.l <- sum(out$time.off > 0)
            #norm.att.off.sq <- (att.off/se.att.off)^2
            #T0.off.p <- 1 - pchisq(sum(norm.att.off.sq[(length(out$time.off) - T0.off.l + 1):length(out$time.off)]), df = T0.off.l)

            att.off.bound <- t(apply(att.off.boot, 1, function(vec) quantile(vec,c(alpha, 1 - alpha), na.rm=TRUE)))
            colnames(att.off.bound) <- c("CI.lower", "CI.upper")
            rownames(att.off.bound) <- out$time.off
        }

        if (!is.null(T.on.carry)) {
            se.carry.att <- apply(carry.att.boot, 1, function(vec) sd(vec, na.rm=TRUE))
            CI.carry.att <- cbind(carry.att - se.carry.att * qnorm(1-alpha/2), 
                                  carry.att + se.carry.att * qnorm(1-alpha/2)) # normal approximation
            pvalue.carry.att <- (1-pnorm(abs(carry.att/se.carry.att)))*2
            est.carry.att <- cbind(carry.att, se.carry.att, 
                             CI.carry.att, pvalue.carry.att)

            colnames(est.carry.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                                      "p.value")
            rownames(est.carry.att) <- carry.time
        }

        if(!is.null(balance.period)){
            se.balance.att <- apply(balance.att.boot, 1, function(vec) sd(vec, na.rm=TRUE))
            CI.balance.att <- cbind(balance.att - se.balance.att * qnorm(1-alpha/2), 
                                    balance.att + se.balance.att * qnorm(1-alpha/2))
            pvalue.balance.att <- (1-pnorm(abs(balance.att/se.balance.att)))*2

            est.balance.att <- cbind(balance.att, se.balance.att, CI.balance.att, 
                                     pvalue.balance.att, out$balance.count)
            colnames(est.balance.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                                           "p.value", "count")
            rownames(est.balance.att) <- out$balance.time
            
            se.balance.avg.att <- sd(balance.avg.att.boot, na.rm=TRUE)
            CI.balance.avg.att <- c(balance.avg.att - se.balance.avg.att  * qnorm(1-alpha/2), 
                                    balance.avg.att + se.balance.avg.att  * qnorm(1-alpha/2))
            p.balance.avg.att <- (1-pnorm(abs(balance.avg.att/se.balance.avg.att)))*2
            est.balance.avg <- t(as.matrix(c(balance.avg.att, se.balance.avg.att, CI.balance.avg.att, p.balance.avg.att)))
            colnames(est.balance.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")

            balance.att.bound <- t(apply(balance.att.boot, 1, function(vec) quantile(vec,c(alpha, 1 - alpha), na.rm=TRUE)))
            colnames(balance.att.bound) <- c("CI.lower", "CI.upper")
            rownames(balance.att.bound) <- out$balance.time

            if (!is.null(placebo.period) & placeboTest == TRUE) {
                balance.att.placebo <- out$balance.att.placebo        
                balance.se.placebo <- sd(balance.att.placebo.boot, na.rm=TRUE)
                balance.CI.placebo <- c(balance.att.placebo - balance.se.placebo * qnorm(1-alpha/2), 
                                        balance.att.placebo + balance.se.placebo * qnorm(1-alpha/2))
                balance.CI.placebo.bound <- c(balance.att.placebo - balance.se.placebo * qnorm(1-alpha), 
                                              balance.att.placebo + balance.se.placebo * qnorm(1-alpha))
                balance.pvalue.placebo <- (1-pnorm(abs(balance.att.placebo/balance.se.placebo)))*2
                est.balance.placebo <- t(as.matrix(c(balance.att.placebo, 
                                            balance.se.placebo, 
                                            balance.CI.placebo, 
                                            balance.pvalue.placebo,
                                            balance.CI.placebo.bound)))
                colnames(est.balance.placebo) <- c("ATT.placebo", "S.E.", 
                                        "CI.lower", "CI.upper", "p.value",
                                        "CI.lower(90%)", "CI.upper(90%)")
            }
        }

        ## average (over time) ATT
        se.avg <- sd(att.avg.boot, na.rm=TRUE)
        CI.avg <- c(att.avg - se.avg * qnorm(1-alpha/2), att.avg + se.avg * qnorm(1-alpha/2))
        pvalue.avg <- (1-pnorm(abs(att.avg/se.avg)))*2
        est.avg <- t(as.matrix(c(att.avg, se.avg, CI.avg, pvalue.avg)))
        colnames(est.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")

        se.avg.unit <- sd(att.avg.unit.boot, na.rm=TRUE)
        CI.avg.unit <- c(att.avg.unit - se.avg.unit * qnorm(1-alpha/2), 
                        att.avg.unit + se.avg.unit * qnorm(1-alpha/2))
        pvalue.avg.unit <- (1-pnorm(abs(att.avg.unit/se.avg.unit)))*2
        est.avg.unit <- t(as.matrix(c(att.avg.unit, se.avg.unit, CI.avg.unit, pvalue.avg.unit)))
        colnames(est.avg.unit) <- c("ATT.avg.unit", "S.E.", "CI.lower", "CI.upper", "p.value")
        
        se.eff.calendar <- apply(calendar.eff.boot, 1, function(vec) sd(vec, na.rm=TRUE))
        CI.eff.calendar <- cbind(calendar.eff - se.eff.calendar * qnorm(1-alpha/2), calendar.eff + se.eff.calendar * qnorm(1-alpha/2))
        pvalue.eff.calendar <- (1-pnorm(abs(calendar.eff/se.eff.calendar)))*2
        est.eff.calendar <- cbind(calendar.eff, se.eff.calendar, CI.eff.calendar, pvalue.eff.calendar,calendar.N)
        colnames(est.eff.calendar) <- c("ATT-calendar", "S.E.", "CI.lower", "CI.upper","p.value", "count")

        se.eff.calendar.fit <- apply(calendar.eff.fit.boot, 1, function(vec) sd(vec, na.rm=TRUE))
        CI.eff.calendar.fit <- cbind(calendar.eff.fit - se.eff.calendar.fit * qnorm(1-alpha/2), calendar.eff.fit + se.eff.calendar.fit * qnorm(1-alpha/2))
        pvalue.eff.calendar.fit <- (1-pnorm(abs(calendar.eff.fit/se.eff.calendar.fit)))*2
        est.eff.calendar.fit <- cbind(calendar.eff.fit, se.eff.calendar.fit, CI.eff.calendar.fit, pvalue.eff.calendar.fit,calendar.N)
        colnames(est.eff.calendar.fit) <- c("ATT-calendar Fitted", "S.E.", "CI.lower", "CI.upper","p.value", "count")

        ## regression coefficents
        if (p > 0) {
            se.beta<-apply(beta.boot, 1, function(vec)sd(vec,na.rm=TRUE))
            CI.beta<-cbind(c(beta) - se.beta * qnorm(1-alpha/2), c(beta) + se.beta * qnorm(1-alpha/2))
            pvalue.beta <- (1-pnorm(abs(beta/se.beta)))*2
            est.beta<-cbind(c(beta), se.beta, CI.beta, pvalue.beta)
            colnames(est.beta)<-c("Coef", "S.E.", "CI.lower", "CI.upper", "p.value")

            if (binary == TRUE) {
                out$marginal[na.pos] <- NA
                se.marginal<-apply(marginal.boot, 1, function(vec)sd(vec,na.rm=TRUE))
                CI.marginal<-cbind(c(out$marginal) - se.marginal * qnorm(1-alpha/2), 
                                c(out$marginal) + se.marginal * qnorm(1-alpha/2))
                pvalue.marginal <- (1-pnorm(abs(out$marginal/se.marginal)))*2
                est.marginal<-cbind(out$marginal, se.marginal, CI.marginal, pvalue.marginal)
                colnames(est.marginal)<-c("marginal", "S.E.", "CI.lower", "CI.upper", "p.value")
            }
        }

        ## placebo test
        if (!is.null(placebo.period) & placeboTest == TRUE) {
            att.placebo <- out$att.placebo        
            se.placebo <- sd(att.placebo.boot, na.rm=TRUE)
            CI.placebo <- c(att.placebo - se.placebo * qnorm(1-alpha/2), 
                            att.placebo + se.placebo * qnorm(1-alpha/2))
            CI.placebo.bound <- c(att.placebo - se.placebo * qnorm(1-alpha), 
                                  att.placebo + se.placebo * qnorm(1-alpha))
            pvalue.placebo <- (1-pnorm(abs(att.placebo/se.placebo)))*2
            est.placebo <- t(as.matrix(c(att.placebo, 
                                         se.placebo, 
                                         CI.placebo, 
                                         pvalue.placebo,
                                         CI.placebo.bound)))
            colnames(est.placebo) <- c("ATT.placebo", "S.E.", 
                                       "CI.lower", "CI.upper", "p.value",
                                       "CI.lower(90%)", "CI.upper(90%)")
        }

        ## carryover test
        if (!is.null(carryover.period) & carryoverTest == TRUE) {
            att.carryover <- out$att.carryover      
            se.carryover <- sd(att.carryover.boot, na.rm=TRUE)
            CI.carryover <- c(att.carryover - se.carryover * qnorm(1-alpha/2), 
                        att.carryover + se.carryover * qnorm(1-alpha/2))
            CI.carryover.bound <- c(att.carryover - se.carryover * qnorm(1-alpha), 
                                    att.carryover + se.carryover * qnorm(1-alpha))
            pvalue.carryover <- (1-pnorm(abs(att.carryover/se.carryover)))*2
            est.carryover <- t(as.matrix(c(att.carryover, se.carryover, 
                                           CI.carryover, pvalue.carryover,
                                           CI.carryover.bound)))
            colnames(est.carryover) <- c("ATT.carryover", "S.E.", 
                                         "CI.lower", "CI.upper", "p.value",
                                         "CI.lower(90%)","CI.upper(90%)")
        }

        ## group effect
        if (!is.null(group)) {
            se.group.att <- apply(group.att.boot, 1, function(vec) sd(vec, na.rm=TRUE))
            CI.group.att <- cbind(c(out$group.att) - se.group.att * qnorm(1-alpha/2), 
                                c(out$group.att) + se.group.att * qnorm(1-alpha/2))
            pvalue.group.att <- (1-pnorm(abs(out$group.att/se.group.att)))*2
            est.group.att <- cbind(out$group.att, se.group.att, CI.group.att, pvalue.group.att)
            colnames(est.group.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value")
        
            est.group.out <- list()
            for(sub.name in group.output.name){
                subgroup.atts <- group.output.origin[[sub.name]]$att.on
                subgroup.atts.boot <- group.atts.boot[[sub.name]]
                subgroup.est.att <- NULL
                subgroup.att.bound <- NULL
                if(dim(subgroup.atts.boot)[1]>0){
                    subgroup.se.att <- apply(subgroup.atts.boot, 1, function(vec) sd(vec, na.rm=TRUE))
                    #subgroup.CI.att <- cbind(subgroup.atts - subgroup.se.att * qnorm(1-alpha/2), 
                    #                        subgroup.atts + subgroup.se.att * qnorm(1-alpha/2)) # normal approximation
                    subgroup.CI.att <- t(apply(subgroup.atts.boot, 1, function(vec) quantile(vec,c(alpha/2, 1 - alpha/2), na.rm=TRUE)))
                    subgroup.pvalue.att <- (1-pnorm(abs(subgroup.atts/subgroup.se.att)))*2
                    subgroup.est.att <- cbind(subgroup.atts, subgroup.se.att , 
                                            subgroup.CI.att, subgroup.pvalue.att, 
                                            group.output.origin[[sub.name]]$count.on)
                    colnames(subgroup.est.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
                                                    "p.value", "count")
                    rownames(subgroup.est.att) <- group.output.origin[[sub.name]]$time.on
                    
                    # for equivalence test
                    subgroup.att.bound <- cbind(subgroup.atts - subgroup.se.att * qnorm(1-alpha), 
                                                subgroup.atts + subgroup.se.att * qnorm(1-alpha)) # one-sided
                    colnames(subgroup.att.bound) <- c("CI.lower", "CI.upper")
                    rownames(subgroup.att.bound) <- group.output.origin[[sub.name]]$time.on
                }
                
                subgroup.att.off.bound <- NULL
                subgroup.est.att.off <- NULL
                if (hasRevs == 1){
                    subgroup.atts.off <- group.output.origin[[sub.name]]$att.off
                    subgroup.atts.off.boot <- group.atts.off.boot[[sub.name]]

                    if(dim(subgroup.atts.off.boot)[1]>0){
                        subgroup.CI.att.off <- t(apply(subgroup.atts.off.boot, 1, function(vec) quantile(vec,c(alpha/2, 1 - alpha/2), na.rm=TRUE)))
                        subgroup.se.att.off <- apply(subgroup.atts.off.boot, 1, function(vec) sd(vec, na.rm=TRUE))
                        subgroup.pvalue.att.off <- apply(subgroup.atts.off.boot, 1, get.pvalue)

                        subgroup.est.att.off <- cbind(subgroup.atts.off, 
                                                    subgroup.se.att.off, 
                                                    subgroup.CI.att.off, 
                                                    subgroup.pvalue.att.off, 
                                                    group.output.origin[[sub.name]]$count.off)
                        colnames(subgroup.est.att.off) <- c("ATT.OFF", "S.E.", "CI.lower", "CI.upper",
                                                            "p.value", "count.off")
                        rownames(subgroup.est.att.off) <- group.output.origin[[sub.name]]$time.off
                        
                        subgroup.att.off.bound <- t(apply(subgroup.atts.off.boot, 1, function(vec) quantile(vec,c(alpha, 1 - alpha), na.rm=TRUE)))
                        colnames(subgroup.att.off.bound) <- c("CI.lower", "CI.upper")
                        rownames(subgroup.att.off.bound) <- group.output.origin[[sub.name]]$time.off                         
                    }   
                }

                ## placebo test
                subgroup.est.placebo <- NULL
                if (!is.null(placebo.period) & placeboTest == TRUE) {
                    subgroup.att.placebo <- group.output.origin[[sub.name]]$att.placebo
                    if(length(subgroup.att.placebo)>0){
                        subgroup.se.placebo <- sd(group.att.placebo.boot[[sub.name]], na.rm=TRUE)
                        subgroup.CI.placebo <- c(subgroup.att.placebo - subgroup.se.placebo * qnorm(1-alpha/2), 
                                                subgroup.att.placebo + subgroup.se.placebo * qnorm(1-alpha/2))
                        subgroup.CI.placebo.bound <- c(subgroup.att.placebo - subgroup.se.placebo * qnorm(1-alpha), 
                                                       subgroup.att.placebo + subgroup.se.placebo * qnorm(1-alpha))
                        subgroup.pvalue.placebo <- (1-pnorm(abs(subgroup.att.placebo/subgroup.se.placebo)))*2
                        subgroup.est.placebo <- t(as.matrix(c(subgroup.att.placebo, 
                                                            subgroup.se.placebo, 
                                                            subgroup.CI.placebo, 
                                                            subgroup.pvalue.placebo,
                                                            subgroup.CI.placebo.bound)))
                        colnames(subgroup.est.placebo) <- c("ATT.placebo", "S.E.", 
                                                            "CI.lower", "CI.upper", "p.value",
                                                            "CI.lower(90%)","CI.upper(90%)")                      
                    }        
                }

                ## carryover test
                subgroup.est.carryover <- NULL
                if (!is.null(carryover.period) & carryoverTest == TRUE) {
                    subgroup.att.carryover <- group.output.origin[[sub.name]]$att.carryover        
                    if(length(subgroup.att.carryover)>0){
                        subgroup.se.carryover <- sd(group.att.carryover.boot[[sub.name]], na.rm=TRUE)
                        subgroup.CI.carryover <- c(subgroup.att.carryover - subgroup.se.carryover * qnorm(1-alpha/2), 
                                                   subgroup.att.carryover + subgroup.se.carryover * qnorm(1-alpha/2))
                        subgroup.CI.carryover.bound <- c(subgroup.att.carryover - subgroup.se.carryover * qnorm(1-alpha), 
                                                         subgroup.att.carryover + subgroup.se.carryover * qnorm(1-alpha))
                        
                        subgroup.pvalue.carryover <- (1-pnorm(abs(subgroup.att.carryover/subgroup.se.carryover)))*2
                        subgroup.est.carryover <- t(as.matrix(c(subgroup.att.carryover, 
                                                            subgroup.se.carryover, 
                                                            subgroup.CI.carryover, 
                                                            subgroup.pvalue.carryover,
                                                            subgroup.CI.carryover.bound)))
                        colnames(subgroup.est.carryover) <- c("ATT.carryover", "S.E.", 
                                                              "CI.lower", "CI.upper", "p.value",
                                                              "CI.lower(90%)","CI.upper(90%)")                 
                    }
                }

                est.group.out[[sub.name]] <- list(att.on = subgroup.est.att,
                                                  att.on.bound = subgroup.att.bound,
                                                  att.on.boot = group.atts.boot[[sub.name]],
                                                  att.off = subgroup.est.att.off,
                                                  att.off.bound = subgroup.att.off.bound,
                                                  att.off.boot = group.atts.off.boot[[sub.name]],
                                                  att.placebo = subgroup.est.placebo,
                                                  att.carryover = subgroup.est.carryover)
            }
        }
    }

    ##storage
    result<-list(est.avg = est.avg,
                 att.avg.boot = att.avg.boot,
                 est.avg.unit = est.avg.unit,
                 att.avg.unit.boot = att.avg.unit.boot,
                 est.eff.calendar = est.eff.calendar,
                 est.eff.calendar.fit = est.eff.calendar.fit,                 
                 est.att = est.att,
                 att.bound = att.bound,
                 att.boot = att.boot,
                 att.count.boot = att.count.boot)

    if (p>0) {
        result <- c(result,list(beta.boot = beta.boot))
        result <- c(result,list(est.beta = est.beta))
        if (binary == TRUE) {
            result <- c(result,list(est.marginal = est.marginal))
        }
    }
    if (hasRevs == 1) {
        result<-c(result,list(est.att.off = est.att.off, 
                              att.off.boot = att.off.boot, 
                              att.off.bound = att.off.bound,
                              att.off.count.boot = att.off.count.boot))
    }

    if (!is.null(T.on.carry)) {
        result <- c(result, list(est.carry.att = est.carry.att))
    }

    if(!is.null(balance.period)){
        result <- c(result, list(est.balance.att = est.balance.att))
        result <- c(result,list(est.balance.avg = est.balance.avg))
        result <- c(result,list(balance.att.bound = balance.att.bound,
                                balance.att.boot = balance.att.boot,
                                balance.count.boot = balance.count.boot))
        if (!is.null(placebo.period) & placeboTest == TRUE) {
            result <- c(result, list(est.balance.placebo = est.balance.placebo, balance.att.placebo.boot = balance.att.placebo.boot))
        }
    } 

    if (!is.null(placebo.period) & placeboTest == TRUE) {
        result <- c(result, list(est.placebo = est.placebo, att.placebo.boot = att.placebo.boot))
    }

    if (!is.null(carryover.period) & carryoverTest == TRUE) {
        result <- c(result, list(est.carryover = est.carryover, att.carryover.boot = att.carryover.boot))
    }

    if (!is.null(group)) {
        result <- c(result, list(est.group.att = est.group.att,
                                 est.group.output = est.group.out))

    }

    return(c(out,result))

    
} ## end of boot


## jackknife se
jackknifed <- function(x,  ## ols estimates
                       y,
                       alpha) { ## sub-sample ols estimates) 

    p <- length(x)
    N <- dim(y)[2]  ## sample size

    X <- matrix(rep(c(x), N), p, N) * N
    Y <- X - y * (N - 1)

    Yvar <- apply(Y, 1, var, na.rm = TRUE)
    vn <- N - apply(is.na(y), 1, sum) 

    Ysd <- sqrt(Yvar/vn)  ## jackknife se

    CI.l <- Ysd * qnorm(alpha/2) + c(x)
    CI.u <- Ysd * qnorm(1 - alpha/2) + c(x)

    ## wald test
    P <- NULL
    for (i in 1:p) {
        subz <- pnorm(c(x)[i]/Ysd[i])
        P <- c(P, 2 * min(1 - subz, subz))
    }

    ## P <- 2 * min(1 - pnorm(c(x)/Ysd), pnorm(c(x)/Ysd))

    out <- list(se = Ysd, CI.l = CI.l, CI.u = CI.u, P = P)

    return(out)
    
}