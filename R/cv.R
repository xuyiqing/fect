###################################################################
## Cross-validation
###################################################################
fect.cv <- function(Y, # Outcome variable, (T*N) matrix
                    X, # Explanatory variables:  (T*N*p) array
                    D, #  Indicator for treated unit (tr==1) 
                    I,
                    II, 
                    T.on, 
                    T.off = NULL, 
                    method = "ife",
                    criterion = "mspe",  
                    k = 5, # CV time
                    cv.prop = 0.1,
                    cv.treat = TRUE, 
                    cv.nobs = 3,
                    r = 0, # initial number of factors considered if CV==1
                    r.end,
                    nlambda = 10, 
                    lambda = NULL,
                    force, 
                    hasRevs = 1,
                    tol, # tolerance level
                    norm.para = NULL,
                    group.level = NULL,
                    group = NULL
                    ) {  
    
    ##-------------------------------##
    ## Parsing data
    ##-------------------------------##  
    placebo.pos <- na.pos <- NULL

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
    oci <- which(c(II) == 1)
    initialOut <- initialFit(data = data.ini, force = force, oci = oci)
    Y0 <- initialOut$Y0
    beta0 <- initialOut$beta0
    if (p > 0 && sum(is.na(beta0)) > 0) {
        beta0[which(is.na(beta0))] <- 0
    }

    ## ------------- restrictions on candidate hyper parameters ---------- ##
    obs.con <- (sum(II) - r.end * (N + TT) + r.end^2 - p) <= 0
    if (obs.con) {
        while((sum(II) - r.end * (N + TT) + r.end^2 - p) <= 0) {
            r.end <- r.end - 1
        }
    }
    if (r.end >= T0.min) {
        if (method %in% c("both", "ife")) {
            cat("Facotr number should not be greater than ", T0.min-1, "\n", sep = "")
        }
        r.end <- T0.min-1
    } else {
        if (obs.con) {
            if (method %in% c("both", "ife")) {
                cat("Facotr number should not be greater than ", r.end, "\n", sep = "")
            }
        }
    }

        ##-------------------------------##
    ## ----------- Main Algorithm ----------- ##
        ##-------------------------------##
    
    validX <- 1 ## no multi-colinearity 
    CV.out.ife <- CV.out.mc <- NULL
    
    ##----------------------------------------------------##
    ##         Cross-validation of r and lambda           ##
    ##----------------------------------------------------##
    
    r.max <- min(TT, r.end)
    r.cv <- 0 ## initial value

    if (method %in% c("ife", "both") && r.max == 0) {
        
        r.cv <- 0
        cat("Cross validation cannot be performed since available pre-treatment records of treated units are too few. So set r.cv = 0.\n ")
        est.best <- inter_fe_ub(YY, Y0, X, II, beta0, 0, force = force, tol)

    } else {

        r.old <- r ## save the minimal number of factors 
        cat("Cross-validating ...","\n") 

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

        ociCV <- matrix(NA, cv.count, k) ## store indicator
        rmCV <- matrix(NA, rm.count, k) ## removed indicator
        Y0CV <- array(NA, dim = c(TT, N, k)) ## store initial Y0
        if (p > 0) {
           beta0CV <- array(NA, dim = c(p, 1, k)) 
        } else {
            beta0CV <- array(0, dim = c(1, 0, k)) ## store initial beta0
        }
        
        ## cv.id.all <- c()
        for (i in 1:k) {
            cv.n <- 0
            repeat{
                cv.n <- cv.n + 1
                ## cv.id <- cv.sample(II, as.integer(sum(II) - cv.count))
                cv.id <- cv.sample(II, D, rm.count, cv.nobs, cv.treat)
                ## cv.id <- sample(oci, as.integer(sum(II) - cv.count), replace = FALSE)
                II.cv <- II
                II.cv[cv.id] <- 0
                con1 <- sum(apply(II.cv, 1, sum) >= 1) == TT
                con2 <- sum(apply(II.cv, 2, sum) >= 1) == N
                if (con1 & con2) {
                    break
                }
                if (cv.n > 100) {
                    stop("Some units have too few pre-treatment observations. Try to remove them.")
                }
                ## cv.id.all <- c(cv.id.all, list(cv.id))
            }
            rmCV[,i] <- cv.id
            ocicv <- setdiff(oci, cv.id)
            ociCV[,i] <- ocicv

            initialOutCv <- initialFit(data = data.ini, force = force, oci = ocicv)
            Y0CV[,,i] <- initialOutCv$Y0
            
            if (p > 0) {
                beta0cv <- initialOutCv$beta0
                if (sum(is.na(beta0cv)) > 0) {
                    beta0cv[which(is.na(beta0cv))] <- 0
                }
                beta0CV[,,i] <- beta0cv
            }
        }
    
        ##  --------------------------------------------- ##
##  ---------------- cross validation for ife model ------------------  ##
        ##  --------------------------------------------- ##
        
        if (method %in% c("ife", "both")) {
            
            # cat("Interactive fixed effects model...\n")
            
            r.pc <- est.pc.best <- MSPE.best <- MSPE.pc.best <- NULL
            if (criterion == "PC") {
                CV.out.ife <- matrix(NA, (r.max - r.old + 1), 6)
                colnames(CV.out.ife) <- c("r", "sigma2", "IC", "PC", "MSPTATT", "MSE")
            } else {
                CV.out.ife <- matrix(NA, (r.max - r.old + 1), 7)
                colnames(CV.out.ife) <- c("r", "sigma2", "IC", "PC", "MSPE", "MSPTATT", "MSE")
            }
            
            CV.out.ife[,"r"] <- c(r.old:r.max)
            CV.out.ife[,"PC"] <- CV.out.ife[,"MSPE"] <- 1e20
                        
            for (i in 1:dim(CV.out.ife)[1]) { ## cross-validation loop starts 

                ## inter FE based on control, before & after 
                r <- CV.out.ife[i, "r"]  
                ## k <- 5
                
                if (criterion %in% c("mspe", "both")) {
                    SSE <- 0
                    for (ii in 1:k) {
                        II.cv <- II
                        II.cv[rmCV[,ii]] <- 0
                        YY.cv <- YY
                        YY.cv[rmCV[,ii]] <- 0
                        est.cv.fit <- inter_fe_ub(YY.cv, as.matrix(Y0CV[,,ii]), X, II.cv, as.matrix(beta0CV[,,ii]), r, force, tol)$fit
                        SSE <- SSE + sum((YY[rmCV[,ii]]-est.cv.fit[rmCV[,ii]])^2)
                    }
                    MSPE <- SSE/(k*(sum(II) - cv.count))
                }
                

                est.cv <- inter_fe_ub(YY, Y0, X, II, beta0, r, force, tol) ## overall
                sigma2 <- est.cv$sigma2 
                IC <- est.cv$IC
                PC <- est.cv$PC

                eff.v.cv <- c(Y - est.cv$fit)[cv.pos]
                meff <- as.numeric(tapply(eff.v.cv, t.on.cv, mean))
                MSPTATT <- sum(meff^2*count.on.cv)/sum(count.on.cv)

                MSE <- sum(eff.v.cv^2)/length(eff.v.cv)

                if(!is.null(norm.para)) {
                    if (criterion %in% c("mspe", "both")) {
                        MSPE <- MSPE*(norm.para[1]^2)
                    }
                    sigma2 <- sigma2*(norm.para[1]^2)
                    IC <- est.cv$IC - log(est.cv$sigma2) + log(sigma2)
                    PC <- PC*(norm.para[1]^2)
                }

                if (criterion %in% c("mspe", "both")) {
                    if ((min(CV.out.ife[,"MSPE"]) - MSPE) > 0.05*min(CV.out.ife[,"MSPE"])) {
                        ## at least 10% improvement for MPSE
                        MSPE.best <- MSPE
                        est.best <- est.cv  
                        r.cv <- r
                    } else {
                        if (r == r.cv + 1) cat("*")
                    }
                }

                if (PC < min(CV.out.ife[,"PC"])) {
                    if (criterion == "both") {
                        MSPE.pc.best <- MSPE
                    }
                    est.pc.best <- est.cv  
                    r.pc <- r
                }


                if (criterion != "pc") {
                    CV.out.ife[i, 2:7] <- c(sigma2, IC, PC, MSPE, MSPTATT, MSE)
                } else {
                    CV.out.ife[i, 2:6] <- c(sigma2, IC, PC, MSPTATT, MSE)
                }
                

                if (criterion == "pc") {
                    cat("\n r = ",r, "; sigma2 = ",
                        sprintf("%.5f",sigma2), "; IC = ",
                        sprintf("%.5f",IC), "; PC = ",
                        sprintf("%.5f",PC), "; MSPTATT = ",
                        sprintf("%.5f",MSPTATT), "; MSE = ",
                        sprintf("%.5f",MSE), sep="")

                } else {
                    cat("\n r = ",r, "; sigma2 = ",
                        sprintf("%.5f",sigma2), "; IC = ",
                        sprintf("%.5f",IC), "; PC = ",
                        sprintf("%.5f",PC), "; MSPE = ",
                        sprintf("%.5f",MSPE), "; MSPTATT = ",
                        sprintf("%.5f",MSPTATT), "; MSE = ",
                        sprintf("%.5f",MSE), sep="")

                }
                
            
            } ## end of while: search for r_star over  

            #MSPE.best <- min(CV.out[,"MSPE"])
            #PC.best <- min(CV.out[,"PC"])

            ## compare 
            if (criterion == "both") {
                if (r.cv > r.pc) {
                    cat("\n\n Factor number selected via cross validation may be larger than the true number. Using the PC criterion.\n\n ")
                    r.cv <- r.pc 
                    est.best <- est.pc.best
                    MSPE.best <- MSPE.pc.best
                }
                est.best.ife <- est.best
                MSPE.best.ife <- MSPE.best
            }
            else if (criterion == "pc") {
                est.best.ife <- est.pc.best  
                r.cv <- r.pc
            }
            else {
                est.best.ife <- est.best
                MSPE.best.ife <- MSPE.best
            }
            
            if (r > (TT-1)) {cat(" (r hits maximum)")}
            cat("\n\n r* = ",r.cv, sep="")
            cat("\n\n") 
        }

        ##  ------------------------------------- ##
##  ---------------- cross validation for mc ---------------  ##
        ##  ------------------------------------- ##
        if (method %in% c("mc", "both")) {
            cat("Matrix completion method...\n")
            eigen.all <- NULL
            if (is.null(lambda) || length(lambda) == 1) {
                ## create the hyper-parameter sequence
                ## biggest candidate lambda 
                ## Y.lambda <- YY 
                ## if (p > 0) {
                ##     for (i in 1:p) {
                ##         Y.lambda <- Y.lambda - X[,,i] * beta0[i,1]
                ##     }
                ## }
                Y.lambda <- YY - Y0
                ## Y.lambda[which(II == 0)] <- Y0[which(II == 0)]
                Y.lambda[which(II == 0)] <- 0
                eigen.all <- svd( Y.lambda / (TT * N) )$d
                lambda.max <- log10(max(eigen.all))
                lambda <- rep(NA, nlambda)
                lambda.by <- 3/(nlambda - 2)
                for (i in 1:(nlambda - 1)) {
                    lambda[i] <- 10^(lambda.max - (i - 1) * lambda.by)
                }
                lambda[nlambda] <- 0
            } else {
                Y.lambda <- YY - Y0
                Y.lambda[which(II == 0)] <- 0
                eigen.all <- svd( Y.lambda / (TT * N) )$d
            }

            ## store all MSPE
            MSPE.best <- NULL 
            CV.out.mc <- matrix(NA, length(lambda), 4)
            colnames(CV.out.mc) <- c("lambda.norm", "MSPE", "MSPTATT", "MSE")
            CV.out.mc[,"lambda.norm"] <- c(lambda/max(eigen.all))
            CV.out.mc[,"MSPE"] <- 1e20

            for (i in 1:length(lambda)) {    
                ## k <- 5
                SSE <- 0
                for (ii in 1:k) {
                    II.cv <- II
                    II.cv[rmCV[,ii]] <- 0
                    YY.cv <- YY
                    YY.cv[rmCV[,ii]] <- 0
                    est.cv.fit <- inter_fe_mc(YY.cv, as.matrix(Y0CV[,,ii]), X, II.cv, as.matrix(beta0CV[,,ii]), 1, lambda[i], force, tol)$fit
                    SSE <- SSE + sum((YY[rmCV[,ii]]-est.cv.fit[rmCV[,ii]])^2)
                }
                MSPE <- SSE/(k*(sum(II) - cv.count))

                est.cv <- inter_fe_mc(YY, Y0, X, II, beta0, 1, lambda[i], force, tol) ## overall
                ## sigma2 <- est.cv$sigma2

                eff.v.cv <- c(Y - est.cv$fit)[cv.pos]
                meff <- as.numeric(tapply(eff.v.cv, t.on.cv, mean))
                MSPTATT <- sum(meff^2*count.on.cv)/sum(count.on.cv) 

                MSE <- sum(eff.v.cv^2)/length(eff.v.cv)

                if(!is.null(norm.para)){
                    MSPE <- MSPE*(norm.para[1]^2)
                    ## sigma2 <- sigma2*(norm.para[1]^2)
                }

                if ((min(CV.out.mc[,"MSPE"]) - MSPE) > 0.05*min(CV.out.mc[,"MSPE"])) {
                    ## at least 10% improvement for MPSE
                    est.best <- est.cv  
                    lambda.cv <- lambda[i]
                    MSPE.best <- MSPE
                } else {
                    if (i > 1) {
                        if (lambda.cv == lambda[i-1]) cat("*")
                    }
                }
                ## CV.out[i, "MSPE"] <- MSPE
                ## CV.out[i, "sigma2"] <- sigma2 
                CV.out.mc[i, 2:4] <- c(MSPE, MSPTATT, MSE)

                cat("\n lambda.norm = ",
                sprintf("%.5f",lambda[i]/max(eigen.all)),"; MSPE = ",
                sprintf("%.5f",MSPE), "; MSPTATT = ",
                sprintf("%.5f",MSPTATT), "; MSE = ", 
                sprintf("%.5f",MSE), sep="")

            }
            est.best.mc <- est.best 
            MSPE.best.mc <- MSPE.best
            cat("\n\n lambda.norm* = ",lambda.cv/max(eigen.all), sep="")
            cat("\n\n")
        }
    }    ## End of Cross-Validation 


    if (method == "ife") {
        est.best <- est.best.ife
        validF <- ifelse(r.cv > 0, 1, 0)
    } 
    else if (method == "mc") {
        est.best <- est.best.mc
        validF <- est.best$validF
    }
    else {
        if (MSPE.best.ife <= MSPE.best.mc) {
            est.best <- est.best.ife
            validF <- ifelse(r.cv > 0, 1, 0) 
            method <- "ife"
        } else {
            est.best <- est.best.mc
            validF <- est.best$validF
            method <- "mc"
        }
        cat("\n\n Recommended method through cross-validation: ", method, sep = "")
        cat("\n\n")
    }    
    validX <- est.best$validX 

        ##------------------------------##
    ## ----------- Summarize -------------- ##
        ##------------------------------## 

    ## 00. run a fect to obtain residuals
    if (method == "ife") {
        if (r.cv == 0) {
            est.fect <- est.best
        } else {
            est.fect <- inter_fe_ub(YY, Y0, X, II, beta0, 0, force = force, tol)
        }
    } else {
        est.fect <- inter_fe_ub(YY, Y0, X, II, beta0, 0, force = force, tol)
    }
           

    ##-------------------------------##
    ##   ATT and Counterfactuals     ##
    ##-------------------------------##

    ## we first adjustment for normalization 
    if (!is.null(norm.para)) {

        Y <- Y * norm.para[1]
        
        if (method == "ife") {
            ## variance of the error term 
            sigma2 <- est.best$sigma2 * (norm.para[1]^2)
            IC <- est.best$IC - log(est.best$sigma2) + log(sigma2)
            PC <- est.best$PC * (norm.para[1]^2)
            est.best$sigma2 <- sigma2
            est.best$IC <- IC
            est.best$PC <- PC
        }

        ## output of estimates
        est.best$mu <- est.best$mu * norm.para[1] 
        if (method == "ife" && r.cv > 0) {
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
        est.fect$fit <- est.fect$fit * norm.para[1]
        est.fect$sigma2 <- est.fect$sigma2 * norm.para[1]
    }

    ## 0. revelant parameters
    sigma2 <- IC <- PC <- NULL
    if (method == "ife") {
        sigma2 <- est.best$sigma2   
        IC <- est.best$IC
        PC <- est.best$PC
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
    eff <- Y - est.best$fit  
    att.avg <- sum(eff * D)/(sum(D))

    ## att.avg.unit
    tr.pos <- which(apply(D, 2, sum) > 0)
    att.unit <- sapply(1:length(tr.pos), function(vec){return(sum(eff[, tr.pos[vec]] * D[, tr.pos[vec]]) / sum(D[, tr.pos[vec]]))})
    att.avg.unit <- mean(att.unit)

    eff.equiv <- Y - est.fect$fit
    equiv.att.avg <- sum(eff.equiv * D) / (sum(D))

    ## 2. rmse for treated units' observations under control
    tr <- which(apply(D, 2, sum) > 0)
    tr.co <- which((as.matrix(1 - D[,tr]) * as.matrix(II[,tr])) == 1)
    eff.tr <- as.matrix(eff[,tr])
    v.eff.tr <- eff.tr[tr.co]
    rmse <- sqrt(mean(v.eff.tr^2))

    ## 3. unbalanced output
    if (0%in%I) {
        eff[which(I == 0)] <- NA
        est.best$fit[which(I == 0)] <- NA
        ## eff.equiv[which(I == 0)] <- NA
    }
    est.best$residuals[which(II == 0)] <- NA
    if (method == "mc") {
        est.best$sigma2 <- mean(c(est.best$residuals[which(II == 1)])^2) ## mean squared error of residuals    
    }    

    ## 4. dynamic effects
    t.on <- c(T.on)
    eff.v <- c(eff) ## a vector
    rm.pos1 <- which(is.na(eff.v))
    rm.pos2 <- which(is.na(t.on))
    eff.v.use1 <- eff.v
    t.on.use <- t.on
    n.on.use <- rep(1:N, each = TT)

    eff.equiv.v <- c(eff.equiv)

    if (NA %in% eff.v | NA %in% t.on) {
        eff.v.use1 <- eff.v[-c(rm.pos1, rm.pos2)]
        t.on.use <- t.on[-c(rm.pos1, rm.pos2)]
        n.on.use <- n.on.use[-c(rm.pos1, rm.pos2)]
        eff.equiv.v <- eff.equiv.v[-c(rm.pos1, rm.pos2)]
    }

    pre.pos <- which(t.on.use <= 0)
    eff.pre <- cbind(eff.v.use1[pre.pos], t.on.use[pre.pos], n.on.use[pre.pos])
    colnames(eff.pre) <- c("eff", "period", "unit")

    eff.pre.equiv <- cbind(eff.equiv.v[pre.pos], t.on.use[pre.pos], n.on.use[pre.pos])
    colnames(eff.pre.equiv) <- c("eff.equiv", "period", "unit")

    pre.sd <- tapply(eff.pre.equiv[,1], eff.pre.equiv[,2], sd)
    pre.sd <- cbind(pre.sd, sort(unique(eff.pre.equiv[, 2])), table(eff.pre.equiv[, 2]))
    colnames(pre.sd) <- c("sd", "period", "count")

    time.on <- sort(unique(t.on.use))
    att.on <- as.numeric(tapply(eff.v.use1, t.on.use, mean)) ## NA already removed
    count.on <- as.numeric(table(t.on.use))

    eff.off <- eff.equiv <- off.sd <- NULL
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

        eff.off.equiv <- cbind(eff.equiv.v[off.pos], t.off.use[off.pos], n.on.use[off.pos])
        colnames(eff.off.equiv) <- c("off.equiv", "period", "unit")

        off.sd <- tapply(eff.off.equiv[,1], eff.off.equiv[,2], sd)
        off.sd <- cbind(off.sd, sort(unique(eff.off.equiv[, 2])), table(eff.off.equiv[, 2]))
        colnames(off.sd) <- c("sd", "period", "count")

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
        sigma2 = est.best$sigma2,
        sigma2.fect = est.fect$sigma2,
        T.on = T.on,
        Y.ct = est.best$fit,
        eff = eff,
        att.avg = att.avg,
        att.avg.unit = att.avg.unit,
        ## supporting
        force = force,
        T = TT,
        N = N,
        p = p, 
        est = est.best,
        method = method,
        mu = est.best$mu,
        beta = beta, 
        validX = validX,
        validF = validF,
        niter = est.best$niter, 
        time = time.on,
        att = att.on,
        count = count.on,
        eff.pre = eff.pre,
        eff.pre.equiv = eff.pre.equiv,
        pre.sd = pre.sd,
        rmse = rmse,
        res = est.best$res)

    if (hasRevs == 1) {
        out <- c(out, list(time.off = time.off, 
                           att.off = att.off,
                           count.off = count.off,
                           eff.off = eff.off,
                           eff.off.equiv = eff.off.equiv,
                           off.sd = off.sd))
    }

    if (force == 1) {
        out<-c(out, list(alpha = est.best$alpha))
    } else if (force == 2) {
        out<-c(out,list(xi = est.best$xi))
    } else if (force == 3) {
        out<-c(out,list(alpha = est.best$alpha, xi = est.best$xi))
    }

    if (method == "ife") {
        out <- c(out, list(r.cv = r.cv, IC = IC, PC = PC))
        if (r.cv > 0) {
            out <- c(out, list(factor = as.matrix(est.best$factor),
                               lambda = as.matrix(est.best$lambda))) 
        }
    }

    if (method == "mc") {
        out <- c(out, list(lambda.cv = lambda.cv, lambda.seq = lambda,
                           lambda.norm = lambda.cv / max(eigen.all), 
                           eigen.all = eigen.all))
    }

    ## CV results
    if (!is.null(CV.out.ife)) {
        out <- c(out, list(CV.out.ife = CV.out.ife))
    }

    if (!is.null(CV.out.mc)) {
        out <- c(out, list(CV.out.mc = CV.out.mc))
    }

    if (!is.null(group)) {
        out <- c(out, list(group.att = group.att))
    }
 
    return(out)
} ## cross-validation function ends

