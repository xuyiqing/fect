###################################################################
## Cross-validation
###################################################################
fect_binary_cv <- function(Y, # Outcome variable, (T*N) matrix
                           X, # Explanatory variables:  (T*N*p) array
                           D, #  Indicator for treated unit (tr==1) 
                           I,
                           II, 
                           T.on, 
                           T.off = NULL, 
                           k = 5, # CV time
                           cv.prop = 0.1,
                           cv.treat = TRUE, 
                           cv.nobs = 3,
                           r = 0, # initial number of factors considered if CV==1
                           r.end,
                           QR = FALSE, 
                           force, 
                           hasRevs = 1,
                           tol, # tolerance level
                           group.level = NULL,
                           group = NULL
                           ) {  
    
    ##-------------------------------##
    ## Parsing data
    ##-------------------------------##  

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
    t.on <- c(T.on)
    T0.min <- min(apply(II, 2, sum))

    ## ----------------- candidate r ------------ ##
    r.max <- min(TT, r.end)
    r.cv <- 0 ## initial value

    ##  --------- initial fit using fastplm --------- ##
    data.ini <- matrix(NA, (TT*N), (2 + 1 + p))
    data.ini[, 2] <- rep(1:N, each = TT)         ## unit fe
    data.ini[, 3] <- rep(1:TT, N)                ## time fe
    data.ini[, 1] <- c(Y)                       ## outcome
    if (p > 0) {                                ## covar
        for (i in 1:p) {
            data.ini[, (3 + i)] <- c(X[, , i])
        }
    }
    ## observed Y0 indicator:
    initialOut <- Y0 <- beta0 <- FE0 <- xi0 <- factor0 <- NULL
    oci <- which(c(II) == 1)


    if (r.max == 0) {
        r.cv <- 0
        message("Cross validation cannot be performed since available pre-treatment records of treated units are too few. So set r.cv = 0.\n ")
        initialOut <- BiInitialFit(data = data.ini, QR = QR, r = 0, force = force, oci = oci)
        Y0 <- initialOut$Y0
        FE0 <- initialOut$FE0
        if (QR == 1) {
            xi0 <- initialOut$xi0
            factor0 <- initialOut$factor0
        }
        if (QR == FALSE) {
            est.best <- inter_fe_d_ub(YY, Y0, FE0, X, II, r.cv, force, tol = tol)
        } else {
            est.best <- inter_fe_d_qr_ub(YY, Y0, FE0, factor0, xi0, X, II, r.cv, force, tol = tol)
        }
    }
    

    ## ------------- restrictions on candidate hyper parameters ---------- ##
    obs.con <- (sum(II) - r.end * (N + TT) + r.end^2 - p) <= 0
    if (obs.con) {
        while((sum(II) - r.end * (N + TT) + r.end^2 - p) <= 0) {
            r.end <- r.end - 1
        }
    }
    if (r.end >= T0.min) {
        message("Facotr number should not be greater than ", T0.min - 1, "\n", sep = "")
        r.end <- T0.min-1
    } else {
        if (obs.con) {
            message("Facotr number should not be greater than ", r.end, "\n", sep = "")
        }
    }

        ##-------------------------------##
    ## ----------- Main Algorithm ----------- ##
        ##-------------------------------##
    
    validX <- 1 ## no multi-colinearity 
    
    ##----------------------------------------------------##
    ##                 Cross-validation of r              ##
    ##----------------------------------------------------##
    
    

    if (r.max > 0) {

        r.old <- r ## save the minimal number of factors 
        message("Cross-validating ...","\n") 

                         ## ----- ##
        ## ------------- initialize ------------ ##
                         ## ----- ##
        
        cv.pos <- which(t.on<=0)
        t.on.cv <- t.on[cv.pos]
        count.on.cv <- as.numeric(table(t.on.cv))
        ## tot.id <- which(c(II)==1) ## observed control data
        ## cv.count <- ceiling((sum(II)*sum(II))/sum(I))
        rm.count <- floor(sum(II)*cv.prop)
        cv.count <- sum(II) - rm.count

        #ociCV <- matrix(NA, cv.count, k) ## store indicator
        #rmCV <- matrix(NA, (length(oci) - cv.count), k) ## removed indicator
        ociCV <- list()
        rmCV <- list()
        Y0CV <- array(NA, dim = c(TT, N, k)) ## store initial Y0
        if (p > 0) {
           beta0CV <- array(NA, dim = c(p, 1, k)) 
        } else {
            beta0CV <- array(0, dim = c(1, 0, k)) ## store initial beta0
        }
        
        for (i in 1:k) {
            cv.n <- 0
            repeat{
                cv.n <- cv.n + 1
                cv.id <- cv.sample(II, D, rm.count, cv.nobs, cv.treat)
                ## cv.id <- cv.sample(II, as.integer(sum(II) - cv.count))
                ## cv.id <- sample(oci, as.integer(sum(II) - cv.count), replace = FALSE)
                #II.cv <- II
                #II.cv[cv.id] <- 0
                II.cv.valid <- II.cv <- II
                II.cv[cv.id] <- 0
                II.cv.valid[cv.id] <- -1
                con1 <- sum(apply(II.cv, 1, sum) >= 1) == TT
                con2 <- sum(apply(II.cv, 2, sum) >= 1) == N
                if (con1 & con2) {
                    break
                }
                if (cv.n>=100) {
                    message("Some units have too few pre-treatment observations. Remove them automatically.")
                    keep.1 <- which(apply(II.cv, 1, sum) < 1)
                    keep.2 <- which(apply(II.cv, 2, sum) < 1)
                    II.cv[keep.1,] <- II[keep.1,]
                    II.cv[,keep.2] <- II[,keep.2]
                    II.cv.valid[keep.1,] <- II[keep.1,]
                    II.cv.valid[,keep.2] <- II[,keep.2]
                    cv.id <- which(II.cv.valid!=II)
                    break
                }
            }
            rmCV[[i]] <- cv.id
            ociCV[[i]] <- setdiff(oci, cv.id)

        }
    
        ##  --------------------------------------------- ##
##  ---------------- cross validation for ife model ------------------  ##
        ##  --------------------------------------------- ##
        
            
        message("Probit model with interactive fixed effects ...\n")
        
        MSPE.best <- NULL
        CV.out <- matrix(NA, (r.max - r.old + 1), 4)
        colnames(CV.out) <- c("r", "IC", "Log-likelihood", "MSPE")
        CV.out[,"r"] <- c(r.old:r.max)
        CV.out[,"MSPE"] <- 1e20

        CVinitialOut <- fit.cv <- Y0.cv <- FE0.cv <- xi0.cv <- factor0.cv <- NULL
                    
        for (i in 1:dim(CV.out)[1]) { ## cross-validation loop starts 

            ## inter FE based on control, before & after 
            r <- CV.out[i, "r"]  
            ## k <- 5
            SSE <- 0
            for (ii in 1:k) {
                II.cv <- II
                II.cv[rmCV[[ii]]] <- 0
                YY.cv <- YY
                YY.cv[rmCV[[ii]]] <- 0 

                CVinitialOut <- BiInitialFit(data = data.ini, QR = QR, r = r, force = force, oci = ociCV[[ii]])
                Y0.cv <- CVinitialOut$Y0
                FE0.cv <- CVinitialOut$FE0
                if (QR == 1) {
                    xi0.cv <- CVinitialOut$xi0
                    factor0.cv <- CVinitialOut$factor0
                }
                if (QR == FALSE) {
                    est.cv.fit <- inter_fe_d_ub(YY.cv, Y0.cv, FE0.cv, X, II.cv, r = r, force, tol = tol)
                } else {
                    est.cv.fit <- inter_fe_d_qr_ub(YY.cv, Y0.cv, FE0.cv, factor0.cv, xi0.cv, X, II.cv, r = r, force, tol = tol)
                }
                fit.cv <- ifelse(est.cv.fit$fit >= 0, 1, 0)

                SSE <- SSE + sum((YY[rmCV[[ii]]]-fit.cv[rmCV[[ii]]])^2)
            }
            MSPE <- SSE/(k*(sum(II) - cv.count))

            ## over-all 
            initialOut <- BiInitialFit(data = data.ini, QR = QR, r = r, force = force, oci = oci)
            Y0 <- initialOut$Y0
            FE0 <- initialOut$FE0
            if (QR == 1) {
                xi0 <- initialOut$xi0
                factor0 <- initialOut$factor0
            }
            if (QR == FALSE) {
                est.cv <- inter_fe_d_ub(YY, Y0, FE0, X, II, r, force, tol = tol)
            } else {
                est.cv <- inter_fe_d_qr_ub(YY, Y0, FE0, factor0, xi0, X, II, r, force, tol = tol)
            }

            IC <- est.cv$IC
            Loglikelihood <- est.cv$loglikelihood

            if (min(CV.out[,"MSPE"]) > MSPE) {
                ## at least 10% improvement for MPSE
                MSPE.best <- MSPE
                est.best <- est.cv  
                r.cv <- r
            } else {
                if (r == r.cv + 1) message("*")
            }

            CV.out[i, 2:4] <- c(IC, Loglikelihood, MSPE)

            message("\n r = ",r, "; IC = ",
                sprintf("%.5f",IC), "; Log-likelihood = ",
                sprintf("%.5f",Loglikelihood), "; MSPE = ",
                sprintf("%.5f",MSPE), sep="")
        }
        if (r > (TT-1)) {message(" (r hits maximum)")}
        message("\n\n r* = ",r.cv, sep="")
        message("\n\n") 

    }    ## End of Cross-Validation 
    validF <- ifelse(r.cv > 0, 1, 0)
    validX <- est.best$validX 

        ##------------------------------##
    ## ----------- Summarize -------------- ##
        ##------------------------------##    

    ##-------------------------------##
    ##   ATT and Counterfactuals     ##
    ##-------------------------------##

    ## 0. relevant parameters
    IC <- est.best$IC
    loglikelihood <- est.best$loglikelihood

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
    Y.ct <- pnorm(est.best$fit)
    eff <- Y - Y.ct    
    att.avg <- sum(eff * D)/(sum(D))

    ## att.avg.unit
    tr.pos <- which(apply(D, 2, sum) > 0)
    att.unit <- sapply(1:length(tr.pos), function(vec){return(sum(eff[, tr.pos[vec]] * D[, tr.pos[vec]]) / sum(D[, tr.pos[vec]]))})
    att.avg.unit <- mean(att.unit)

    ## 2. rmse for treated units' observations under control: useless for binary
    ## average marginal effect
    marginal <- NULL
    #if (binary == TRUE) {
        if (p > 0) {
            dense <- dnorm(c(est.best$fit[which(II == 1)]))
            marginal <- as.matrix(sapply(1:p, function(vec){ mean(beta[vec] * dense)} ))
        }
    #}

    ## 3. unbalanced output
    if (0%in%I) {
        eff[which(I == 0)] <- NA
        Y.ct[which(I == 0)] <- NA
        est.best$fit[which(I == 0)] <- NA
    } 

    ## 4. dynamic effects
    t.on <- c(T.on)
    eff.v <- c(eff) ## a vector
    rm.pos1 <- which(is.na(eff.v))
    rm.pos2 <- which(is.na(t.on))
    eff.v.use1 <- eff.v
    t.on.use <- t.on
    n.on.use <- rep(1:N, each = TT)

    if (NA %in% eff.v | NA %in% t.on) {
        eff.v.use1 <- eff.v[-c(rm.pos1, rm.pos2)]
        t.on.use <- t.on[-c(rm.pos1, rm.pos2)]
        n.on.use <- n.on.use[-c(rm.pos1, rm.pos2)]
    }

    pre.pos <- which(t.on.use <= 0)
    eff.pre <- cbind(eff.v.use1[pre.pos], t.on.use[pre.pos], n.on.use[pre.pos])
    colnames(eff.pre) <- c("eff", "period", "unit")

    time.on <- sort(unique(t.on.use))
    att.on <- as.numeric(tapply(eff.v.use1, t.on.use, mean)) ## NA already removed
    count.on <- as.numeric(table(t.on.use))

    ## 5. placebo effect, if placeboTest == 1 , no placebo for CV


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

        time.off <- sort(unique(t.off.use))
        att.off <- as.numeric(tapply(eff.v.use2, t.off.use, mean)) ## NA already removed
        count.off <- as.numeric(table(t.off.use))
    }

    ## 7. cohort effects
    if (!is.null(group)) {
        cohort <- cbind(c(group), c(D), c(eff.v))
        rm.pos <- unique(c(rm.pos1, which(cohort[, 2] == 0)))
        cohort <- cohort[-rm.pos, ]

        g.level <- sort(unique(cohort[, 1]))
        raw.group.att <- as.numeric(tapply(cohort[, 3], cohort[, 1], mean))

        group.att <- rep(NA, length(group.level))
        group.att[which(group.level %in% g.level)] <- raw.group.att 
    }

    ##-------------------------------##
    ##           Storage 
    ##-------------------------------##  

    ##control group residuals
    out<-list(
        ## main results 
        T.on = T.on,
        Y.ct = Y.ct,
        eff = eff,
        att.avg = att.avg,
        att.avg.unit = att.avg.unit,
        ## supporting
        force = force,
        T = TT,
        N = N,
        p = p, 
        r.cv = r.cv,
        est = est.best,
        mu = est.best$mu,
        beta = beta, 
        validX = validX,
        validF = validF,
        niter = est.best$niter, 
        time = time.on,
        att = att.on,
        count = count.on,
        eff.pre = eff.pre,
        CV.out = CV.out
    )

    if (hasRevs == 1) {
        out <- c(out, list(time.off = time.off, 
                           att.off = att.off,
                           count.off = count.off))
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

    if (!is.null(group)) {
        out <- c(out, list(group.att = group.att))
    }
 
    return(out)
} ## binary cross-validation function ends

