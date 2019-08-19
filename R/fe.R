###################################################################
## IFE Model Function
###################################################################
fect.fe <- function(Y, # Outcome variable, (T*N) matrix
                    X, # Explanatory variables:  (T*N*p) array
                    D, #  Indicator for treated unit (tr==1) 
                    I,
                    II, 
                    T.on, 
                    T.off = NULL, 
                    r.cv = 0, # initial number of factors considered if CV==1
                    binary = FALSE,
                    QR = FALSE, 
                    force, 
                    hasRevs = 1,
                    tol, # tolerance level
                    boot = FALSE, # bootstrapped sample
                    placeboTest = 0,
                    placebo.period = NULL,
                    norm.para = NULL,
                    time.on.seq = NULL,
                    time.off.seq = NULL) {  
    
    ##-------------------------------##
    ## Parsing data
    ##-------------------------------##  
    placebo.pos <- na.pos <- NULL
    res.sd1 <- res.sd2 <- NULL

    ## unit id and time
    TT <- dim(Y)[1]
    N <- dim(Y)[2]
    if (is.null(X) == FALSE) {
        p <- dim(X)[3]
    } else {
        p <- 0
        X <- array(0, dim = c(1, 1, 0))
    }

    ## replicate data
    YY <- Y
    YY[which(II == 0)] <- 0 ## reset to 0 

    ## initial fit using fastplm
    data.ini <- matrix(NA, (TT*N), (2 + 1 + p))
    data.ini[, 2] <- rep(1:N, each = TT)         ## unit fe
    data.ini[, 3] <- rep(1:TT, N)                ## time fe
    data.ini[, 1] <- c(Y)                        ## outcome
    if (p > 0) {                                 ## covar
        for (i in 1:p) {
            data.ini[, (3 + i)] <- c(X[, , i])
        }
    }
    ## observed Y0 indicator:
    initialOut <- Y0 <- beta0 <- FE0 <- xi0 <- factor0 <- NULL
    oci <- which(c(II) == 1)
    if (binary == FALSE) { 
        initialOut <- initialFit(data = data.ini, force = force, oci = oci)
        Y0 <- initialOut$Y0
        beta0 <- initialOut$beta0
        if (p > 0 && sum(is.na(beta0)) > 0) {
            beta0[which(is.na(beta0))] <- 0
        }
        ## ini.res <- initialOut$res
    } else {
        initialOut <- BiInitialFit(data = data.ini, QR = QR, r = r.cv, force = force, oci = oci)
        Y0 <- initialOut$Y0
        beta0 <- initialOut$beta0
        FE0 <- initialOut$FE0
        if (QR == 1) {
            xi0 <- initialOut$xi0
            factor0 <- initialOut$factor0
        }
    }
 
        ##-------------------------------##
    ## ----------- Main Algorithm ----------- ##
        ##-------------------------------##
    
    validX <- 1 ## no multi-colinearity
    est.equiv <- NULL
    if (binary == FALSE) {
        est.best <- inter_fe_ub(YY, Y0, X, II, beta0, r.cv, force = force, tol) 
        if (boot == FALSE) {
            if (r.cv == 0) {
                est.equiv <- est.best
            } else {
                est.equiv <- inter_fe_ub(YY, Y0, X, II, beta0, 0, force = force, tol)
            }
        }
    } else {
        if (QR == FALSE) {
            est.best <- inter_fe_d_ub(YY, Y0, FE0, X, II, r.cv, force, tol = tol)
        } else {
            est.best <- inter_fe_d_qr_ub(YY, Y0, FE0, factor0, xi0, X, II, r.cv, force, tol = tol)
        }
        
    }
    validX <- est.best$validX
    validF <- ifelse(r.cv > 0, 1, 0)
        
        ##------------------------------##
    ## ----------- Summarize -------------- ##
        ##------------------------------##    

    ##-------------------------------##
    ##   ATT and Counterfactuals     ##
    ##-------------------------------##

    ## we first adjustment for normalization 
    if (!is.null(norm.para) && binary == FALSE) {
        Y <- Y * norm.para[1]
        ## variance of the error term 
        sigma2 <- est.best$sigma2 * (norm.para[1]^2)
        IC <- est.best$IC - log(est.best$sigma2) + log(sigma2)
        PC <- est.best$PC * (norm.para[1]^2)
        est.best$sigma2 <- sigma2
        est.best$IC <- IC
        est.best$PC <- PC

        ## output of estimates
        est.best$mu <- est.best$mu * norm.para[1] 
        if (r.cv > 0) {
            est.best$lambda <- est.best$lambda * norm.para[1]
            est.best$VNT <- est.best$VNT * norm.para[1]
        }
        if (force%in%c(1, 3)) {
            est.best$alpha <- est.best$alpha * norm.para[1] 
        }
        if (force%in%c(2,3)) {
            est.best$xi <- est.best$xi * norm.para[1] 
        }
        #if (p>0) {
        #    est.best$beta <- est.best$beta * norm.para[1]
        #}
        est.best$residuals <- est.best$residuals * norm.para[1] 
        est.best$fit <- est.best$fit * norm.para[1]
        ## ini.res <- ini.res * norm.para[1] 
        if (boot == FALSE) {
            est.equiv$fit <- est.equiv$fit * norm.para[1]
        }
    }

    ## 0. relevant parameters
    IC <- est.best$IC
    if (binary == FALSE) {
        sigma2 <- est.best$sigma2   
        PC <- est.best$PC
    } else {
        loglikelihood <- est.best$loglikelihood
    }

    if (p>0) {
        na.pos <- is.nan(est.best$beta)
        beta <- est.best$beta
        if( sum(na.pos) > 0 ) {
            beta[na.pos] <- NA
        }
    } else {
        beta <- NA
    }
   
    ## 1. estimated att and counterfactuals
    Y.ct.equiv <- Y.ct <- NULL
    if (binary == FALSE) {
        Y.ct <- est.best$fit
        if (boot == FALSE) {
            Y.ct.equiv <- est.equiv$fit
        }
    } else {
        Y.ct <- pnorm(est.best$fit)
    }
    eff <- Y - Y.ct    
    att.avg <- sum(eff * D)/(sum(D))
    
    equiv.att.avg <- eff.equiv <- NULL
    if (binary == FALSE && boot == FALSE) {
        eff.equiv <- Y - Y.ct.equiv
        equiv.att.avg <- sum(eff.equiv * D)/(sum(D))
    }

    ## 2. rmse for treated units' observations under control
    if (binary == 0) {
        tr <- which(apply(D, 2, sum) > 0)
        tr.co <- which((as.matrix(1 - D[,tr]) * as.matrix(II[,tr])) == 1)
        eff.tr <- as.matrix(eff[,tr])
        v.eff.tr <- eff.tr[tr.co]
        rmse <- sqrt(mean(v.eff.tr^2))
    }

    ## 3. unbalanced output
    if (0 %in% I) {
        eff[which(I == 0)] <- NA
        Y.ct[which(I == 0)] <- NA
        est.best$fit[which(I == 0)] <- NA
    }
    if (binary == FALSE) {
        est.best$residuals[which(II == 0)] <- NA 
    }
      

    ## 4. dynamic effects
    t.on <- c(T.on)
    eff.v <- c(eff) ## a vector
    
    eff.equiv.v <- NULL
    if (binary == FALSE && boot == FALSE) {
        eff.equiv.v <- c(eff.equiv)
    }

    rm.pos1 <- which(is.na(eff.v))
    rm.pos2 <- which(is.na(t.on)) 

    eff.v.use1 <- eff.v
    t.on.use <- t.on
    n.on.use <- rep(1:N, each = TT)

    if (NA %in% eff.v | NA %in% t.on) {
        eff.v.use1 <- eff.v[-c(rm.pos1, rm.pos2)]
        t.on.use <- t.on[-c(rm.pos1, rm.pos2)]
        n.on.use <- n.on.use[-c(rm.pos1, rm.pos2)]
        if (binary == FALSE && boot == FALSE) {
            eff.equiv.v <- eff.equiv.v[-c(rm.pos1, rm.pos2)]
        }
    }

    pre.pos <- which(t.on.use <= 0)
    eff.pre <- cbind(eff.v.use1[pre.pos], t.on.use[pre.pos], n.on.use[pre.pos])
    colnames(eff.pre) <- c("eff", "period", "unit")

    pre.sd <- eff.pre.equiv <- NULL
    if (binary == FALSE && boot == FALSE) {
        eff.pre.equiv <- cbind(eff.equiv.v[pre.pos], t.on.use[pre.pos], n.on.use[pre.pos])
        colnames(eff.pre.equiv) <- c("eff.equiv", "period", "unit")

        pre.sd <- tapply(eff.pre.equiv[,1], eff.pre.equiv[,2], sd)
        pre.sd <- cbind(pre.sd, sort(unique(eff.pre.equiv[, 2])), table(eff.pre.equiv[, 2]))
        colnames(pre.sd) <- c("sd", "period", "count")
    }

    time.on <- sort(unique(t.on.use))
    att.on <- as.numeric(tapply(eff.v.use1, t.on.use, mean)) ## NA already removed
    count.on <- as.numeric(table(t.on.use))

    if (!is.null(time.on.seq)) {
        count.on.med <- att.on.med <- rep(NA, length(time.on.seq))
        att.on.med[which(time.on.seq %in% time.on)] <- att.on
        count.on.med[which(time.on.seq %in% time.on)] <- count.on
        att.on <- att.on.med
        count.on <- count.on.med
        time.on <- time.on.seq
    }

    ## 5. placebo effect, if placeboTest == 1 
    if (!is.null(placebo.period) && placeboTest == 1) {              
        if (length(placebo.period) == 1) {
            placebo.pos <- which(time.on == placebo.period)
            att.placebo <- att.on[placebo.pos]
        } else {
            placebo.pos <- which(time.on >= placebo.period[1] & time.on <= placebo.period[2])
            att.placebo <- sum(att.on[placebo.pos] * count.on[placebo.pos]) / sum(count.on[placebo.pos])
        }
    }

    eff.off.equiv <- off.sd <- eff.off <- NULL

    ## 6. switch-off effects
    if (hasRevs == 1) {    
        t.off <- c(T.off)
        rm.pos3 <- which(is.na(t.off))
        eff.v.use2 <- eff.v
        t.off.use <- t.off

        if (NA %in% eff.v | NA %in% t.off) {
            eff.v.use2 <- eff.v[-c(rm.pos1, rm.pos3)]
            t.off.use <- t.off[-c(rm.pos1, rm.pos3)]
        }

        off.pos <- which(t.off.use > 0)
        eff.off <- cbind(eff.v.use2[off.pos], t.off.use[off.pos], n.on.use[off.pos])
        colnames(eff.off) <- c("eff", "period", "unit")

        if (binary == FALSE && boot == FALSE) {
            eff.off.equiv <- cbind(eff.equiv.v[off.pos], t.off.use[off.pos], n.on.use[off.pos])
            colnames(eff.off.equiv) <- c("off.equiv", "period", "unit")

            off.sd <- tapply(eff.off.equiv[,1], eff.off.equiv[,2], sd)
            off.sd <- cbind(off.sd, sort(unique(eff.off.equiv[, 2])), table(eff.off.equiv[, 2]))
            colnames(off.sd) <- c("sd", "period", "count")
        }

        time.off <- sort(unique(t.off.use))
        att.off <- as.numeric(tapply(eff.v.use2, t.off.use, mean)) ## NA already removed
        count.off <- as.numeric(table(t.off.use))

        if (!is.null(time.off.seq)) {
            count.off.med <- att.off.med <- rep(NA, length(time.off.seq))
            att.off.med[which(time.off.seq %in% time.off)] <- att.off
            count.off.med[which(time.off.seq %in% time.off)] <- count.off
            att.off <- att.off.med
            count.off <- count.off.med
            time.off <- time.off.seq
        }
    }

    method <- ifelse(r.cv > 0, "ife", "fe")

  
    ##-------------------------------##
    ##            Storage            ##
    ##-------------------------------##  
    out<-list(
        ## main results 
        method = method,
        Y.ct = Y.ct,
        eff = eff,
        att.avg = att.avg,
        ## supporting
        force = force,
        T = TT,
        N = N,
        p = p,
        r.cv = r.cv, 
        IC = IC, 
        beta = beta,
        est = est.best,
        mu = est.best$mu,
        niter = est.best$niter,
        validX = validX,
        validF = validF,
        time = time.on,
        att = att.on,
        count = count.on,
        eff.pre = eff.pre,
        eff.pre.equiv = eff.pre.equiv,
        pre.sd = pre.sd)

    if (binary == 0) {
        out <- c(out, list(PC = PC,
                           sigma2 = sigma2, 
                           res = est.best$residuals,
                           rmse = rmse))
        #if (boot == FALSE) {
        #    out <- c(out, list(equiv.att.avg = equiv.att.avg))
        #}
    } else {
        out <- c(out, list(loglikelihood = loglikelihood))
    }
    
    if (hasRevs == 1) {
        out <- c(out, list(time.off = time.off, 
                           att.off = att.off,
                           count.off = count.off,
                           eff.off = eff.off,
                           eff.off.equiv = eff.off.equiv,
                           off.sd = off.sd))
    }
    if (r.cv > 0) {
        out<-c(out,list(factor = as.matrix(est.best$factor),
                        lambda = as.matrix(est.best$lambda))) 
    }

    if (force == 1) {
        out<-c(out, list(alpha = est.best$alpha))
    } else if (force == 2) {
        out<-c(out,list(xi = est.best$xi))
    } else if (force == 3) {
        out<-c(out,list(alpha = est.best$alpha, xi = est.best$xi))
    }

    if (!is.null(placebo.period) && placeboTest == 1) {
        out <- c(out, list(att.placebo = att.placebo))
    }

    return(out)
} ## fe functions ends. 