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
                      method = "fe",
                      criterion = "mspe",
                      CV,
                      k = 5,
                      r = 0, 
                      r.end,
                      lambda = NULL,
                      nlambda = 10,
                      binary,
                      QR,
                      force,                      
                      hasRevs = 1,
                      tol,
                      norm.para,
                      placebo.period = NULL,
                      placeboTest = FALSE,
                      nboots,
                      parallel = TRUE,
                      cores = NULL) {
    
    
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
    } else {
        ## treatement indicator
        tr <- which(apply(D, 2, sum) > 0)
        co <- which(apply(D, 2, sum) == 0)

        Ntr <- length(tr)
        Nco <- length(co)
    }

    
    ## estimation
    if (CV == 0) { 
        if (method == "fe") {
            out <- fect.fe(Y = Y, X = X, D = D, I = I, II = II, 
                           T.on = T.on, T.off = T.off,
                           r.cv = r, binary = binary, QR = QR,
                           force = force, hasRevs = hasRevs, 
                           tol = tol, boot = 0,
                           norm.para = norm.para, 
                           placebo.period = placebo.period,
                           placeboTest = placeboTest)
        } else {
            out <- fect.mc(Y = Y, X = X, D = D, I = I, II = II,
                           T.on = T.on, T.off = T.off, 
                           lambda.cv = lambda, force = force, hasRevs = hasRevs, 
                           tol = tol, boot = 0,
                           norm.para = norm.para,
                           placebo.period = placebo.period,
                           placeboTest = placeboTest)
        }
    } else {
        ## cross-valiadtion 
        if (binary == 0) {
            out <- fect.cv(Y = Y, X = X, D = D, I = I, II = II, 
                       T.on = T.on, T.off = T.off,
                       method = method, criterion = criterion,
                       k = k, r = r, r.end = r.end, 
                       nlambda = nlambda, lambda = lambda, 
                       force = force, hasRevs = hasRevs, 
                       tol = tol, norm.para = norm.para)

            method <- out$method
        } else {
            out <- fect.binary.cv(Y = Y, X = X, D = D, 
                                  I = I, II = II, 
                                  T.on = T.on, T.off = T.off,
                                  k = k, r = r, r.end = r.end, 
                                  QR = QR, force = force, 
                                  hasRevs = hasRevs, tol = tol)
            method <- "fe"
        }
        
    }
    
    
    ## output
    validX <- out$validX
    eff <- out$eff
    att.avg <- out$att.avg

    att.on <- out$att.on
    time.on <- out$time.on

    time.off <- NULL
    if (hasRevs == 1) {
        att.off <- out$att.off
        time.off <- out$time.off
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
 
    ## bootstrapped estimates
    eff.boot <- array(0,dim = c(TT, Ntr, nboots))  ## to store results
    att.avg.boot <- matrix(0, nboots, 1)
    att.on.boot <- matrix(0, length(time.on), nboots)
    att.on.count.boot <- matrix(0, length(time.on), nboots)
    if (hasRevs == 1) {
        att.off.boot <- matrix(0, length(time.off), nboots) 
        att.off.count.boot <- matrix(0, length(time.off), nboots)   
    }
    if (p > 0) {
        beta.boot <- matrix(0, p, nboots)
    }
    if (!is.null(placebo.period) & placeboTest == TRUE) {
        att.placebo.boot <- matrix(0, nboots, 1)
    } 

    cat("\rBootstrapping ...\n")
 
    if (method == "fe") {
        one.nonpara <- function() {

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
                    } else {
                        repeat{
                            boot.id <- sample(tr, Ntr, replace=TRUE)
                            if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
                                break
                            }
                        }
                    }
                } else {
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
                        } else {
                            repeat{
                                fake.tr <- sample(tr, Ntr, replace=TRUE)
                                fake.rev <- sample(rev, Nrev, replace=TRUE)
                                boot.id <- c(fake.rev, fake.tr)
                                if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
                                    break
                                }
                            }
                        }
                    } else {
                        if (Nco > 0) {
                            repeat{
                                fake.co <- sample(co, Nco, replace=TRUE)
                                fake.rev <- sample(rev, Nrev, replace=TRUE)
                                boot.id <- c(fake.rev, fake.co)
                                if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
                                    break
                                }
                            }
                        } else {
                            repeat{
                                boot.id <- sample(rev, Nrev, replace=TRUE)
                                if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
                                    break
                                }
                            }
                        }
                    }
                }
            } else {
                cl.boot <- sample(cl.unique, length(cl.unique), replace = TRUE)
                cl.boot.uni <- unique(cl.boot)
                cl.boot.count <- as.numeric(table(cl.boot))
                boot.id <- c()
                for (kk in 1:length(cl.boot.uni)) {
                    boot.id <- c(boot.id, rep(which(cl == cl.boot.uni[kk]), cl.boot.count[kk]))
                }
            }
            X.boot <- X[,boot.id,,drop = FALSE]
            D.boot <- D[, boot.id]
            I.boot <- I[, boot.id]

            if (sum(c(D.boot) == 0) == 0 | sum(c(D.boot) == 1) == 0 | sum(c(I.boot) == 1) == 0) {
                boot0 <- list(att.avg = NA, att.on = NA, count.on = NA, 
                              beta = NA, att.off = NA, count.off = NA, 
                              att.placebo = NA)
                return(boot0)
            } else {
                T.off.boot <- NULL
                if (hasRevs == TRUE) {
                    T.off.boot <- T.off[, boot.id]
                }
                placebo.period.boot <- NULL
                if (placeboTest == TRUE) {
                    placebo.period.boot <- placebo.period
                }
                    
                boot <- try(fect.fe(Y = Y[, boot.id], X = X.boot, D = D.boot,
                                    I = I.boot, II = II[, boot.id], 
                                    T.on = T.on[, boot.id], T.off = T.off.boot, 
                                    r.cv = out$r.cv, binary = binary,
                                    QR = QR, force = force,
                                    hasRevs = hasRevs, tol = tol, boot = 1,
                                    norm.para = norm.para,
                                    time.on.seq = time.on, time.off.seq = time.off,
                                    placebo.period = placebo.period.boot, 
                                    placeboTest = placeboTest), silent = TRUE)

                if ('try-error' %in% class(boot)) {
                    boot0 <- list(att.avg = NA, att.on = NA, count.on = NA, 
                                  beta = NA, att.off = NA, count.off = NA, 
                                  att.placebo = NA)
                    return(boot0)
                } else {
                    return(boot)
                }
            }            
        } 
            
    } else { ## mc
        one.nonpara <- function() {

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
                    } else {
                        repeat{
                            boot.id <- sample(tr, Ntr, replace=TRUE)
                            if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
                                break
                            }
                        }
                    }
                } else {
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
                        } else {
                            repeat{
                                fake.tr <- sample(tr, Ntr, replace=TRUE)
                                fake.rev <- sample(rev, Nrev, replace=TRUE)
                                boot.id <- c(fake.rev, fake.tr)
                                if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
                                    break
                                }
                            }
                        }
                    } else {
                        if (Nco > 0) {
                            repeat{
                                fake.co <- sample(co, Nco, replace=TRUE)
                                fake.rev <- sample(rev, Nrev, replace=TRUE)
                                boot.id <- c(fake.rev, fake.co)
                                if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
                                    break
                                }
                            }
                        } else {
                            repeat{
                                boot.id <- sample(rev, Nrev, replace=TRUE)
                                if (sum(apply(as.matrix(I[,boot.id]),1,sum)>=1)==TT) {
                                    break
                                }
                            }
                        }
                    }
                }
            } else {
                cl.boot <- sample(cl.unique, length(cl.unique), replace = TRUE)
                cl.boot.uni <- unique(cl.boot)
                cl.boot.count <- as.numeric(table(cl.boot))
                boot.id <- c()
                for (kk in 1:length(cl.boot.uni)) {
                    boot.id <- c(boot.id, rep(which(cl == cl.boot.uni[kk]), cl.boot.count[kk]))
                }
            }
            
                
            X.boot <- X[,boot.id,,drop = FALSE]
            D.boot <- D[, boot.id]
            I.boot <- I[, boot.id]

            if (sum(c(D.boot) == 0) == 0 | sum(c(D.boot) == 1) == 0 | sum(c(I.boot) == 1) == 0) {
                boot0 <- list(att.avg = NA, att.on = NA, count.on = NA, 
                              beta = NA, att.off = NA, count.off = NA, 
                              att.placebo = NA)
                return(boot0)
            } else {
                T.off.boot <- NULL
                if (hasRevs == TRUE) {
                    T.off.boot <- T.off[, boot.id]
                }
                placebo.period.boot <- NULL
                if (placeboTest == TRUE) {
                    placebo.period.boot <- placebo.period
                }
                
                boot <- try(fect.mc(Y = Y[,boot.id], X = X.boot, D = D[,boot.id],
                                    I = I[,boot.id], II = II[,boot.id],
                                    T.on = T.on[,boot.id], T.off = T.off.boot, 
                                    lambda.cv = out$lambda.cv, force = force, 
                                    hasF = out$validF, hasRevs = hasRevs, 
                                    tol = tol, boot = 1,
                                    norm.para = norm.para,
                                    time.on.seq = time.on, time.off.seq = time.off,
                                    placebo.period = placebo.period.boot, 
                                    placeboTest = placeboTest), silent = TRUE)
                
                if ('try-error' %in% class(boot)) {
                    boot0 <- list(att.avg = NA, att.on = NA, count.on = NA, 
                                  beta = NA, att.off = NA, count.off = NA, 
                                  att.placebo = NA)
                    return(boot0)
                } else {
                    return(boot)
                }
            }                        
        } 
    }
    
    ## computing
    if (parallel == TRUE) { 
        boot.out <- foreach(j=1:nboots, 
                            .inorder = FALSE,
                            .export = c("fect.fe", "fect.mc", "get_term"),
                            .packages = c("fect")
                            ) %dopar% {
                                return(one.nonpara())
                            }

        for (j in 1:nboots) { 
            att.avg.boot[j,] <- boot.out[[j]]$att.avg
            att.on.boot[,j] <- boot.out[[j]]$att.on 
            att.on.count.boot[,j] <- boot.out[[j]]$count.on 
            if (p > 0) {
                beta.boot[,j] <- boot.out[[j]]$beta
            }
            if (hasRevs == 1) {
                att.off.boot[,j] <- boot.out[[j]]$att.off
                att.off.count.boot[,j] <- boot.out[[j]]$count.off 
            }
            if (!is.null(placebo.period) & placeboTest == TRUE) {
                att.placebo.boot[j,] <- boot.out[[j]]$att.placebo
            }
        } 
    } else {
        for (j in 1:nboots) { 
            boot <- one.nonpara() 
            att.avg.boot[j,] <- boot$att.avg
            att.on.boot[,j] <- boot$att.on
            att.on.count.boot[,j] <- boot$count.on 
            if (p > 0) {
                beta.boot[,j] <- boot$beta
            }
            if (hasRevs == 1) {
                att.off.boot[,j] <- boot$att.off
                att.off.count.boot[,j] <- boot$count.off 
            }
            if (!is.null(placebo.period) & placeboTest == TRUE) {
                att.placebo.boot[j,] <- boot$att.placebo
            }
            ## report progress
            if (j%%100 == 0)  {
                cat(".")   
            }  
        }  
    } 
    ## end of bootstrapping
    cat("\r")

    ## remove failure bootstrap
    ## alternative condition? max(apply(is.na(att.on.boot),2,sum)) == dim(att.on.boot)[1]
    if (sum(is.na(c(att.avg.boot))) > 0) {
        boot.rm <- which(is.na(c(att.avg.boot)))
        att.avg.boot <- as.matrix(att.avg.boot[-boot.rm,])
        att.on.boot <- as.matrix(att.on.boot[,-boot.rm])
        att.on.count.boot <- as.matrix(att.on.count.boot[,-boot.rm])
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
        if (!is.null(placebo.period) & placeboTest == TRUE) {
            att.placebo.boot <- as.matrix(att.placebo.boot[-boot.rm,])
        }

    }
    cat("Actual bootstrap times: ", dim(att.on.boot)[2], "\n", sep = "")
     
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
    CI.att.on <- t(apply(att.on.boot, 1, function(vec) quantile(vec,c(0.025,0.975), na.rm=TRUE)))
    se.att.on <- apply(att.on.boot, 1, function(vec) sd(vec, na.rm=TRUE))
    pvalue.att.on <- apply(att.on.boot, 1, get.pvalue)

    est.att.on <- cbind(att.on, se.att.on, CI.att.on, pvalue.att.on, out$count.on)
    colnames(est.att.on) <- c("ATT.ON", "S.E.", "CI.lower", "CI.upper",
                              "p.value", "count.on")
    rownames(est.att.on) <- out$time.on
    T0.on.l <- sum(out$time.on <= 0)
    norm.att.on.sq <- (att.on/se.att.on)^2
    T0.on.p <- 1 - pchisq(sum(norm.att.on.sq[1:T0.on.l]), df = T0.on.l)

    if (hasRevs == 1) {
        CI.att.off <- t(apply(att.off.boot, 1, function(vec) quantile(vec,c(0.025,0.975), na.rm=TRUE)))
        se.att.off <- apply(att.off.boot, 1, function(vec) sd(vec, na.rm=TRUE))
        pvalue.att.off <- apply(att.off.boot, 1, get.pvalue)

        est.att.off <- cbind(att.off, se.att.off, CI.att.off, pvalue.att.off, out$count.off)
        colnames(est.att.off) <- c("ATT.OFF", "S.E.", "CI.lower", "CI.upper",
                                   "p.value", "count.off")
        rownames(est.att.off) <- out$time.off
        T0.off.l <- sum(out$time.off > 0)
        norm.att.off.sq <- (att.off/se.att.off)^2
        T0.off.p <- 1 - pchisq(sum(norm.att.off.sq[(length(out$time.off) - T0.off.l + 1):length(out$time.off)]), df = T0.off.l)
    }

    ## average (over time) ATT
    CI.avg <- quantile(att.avg.boot, c(0.025,0.975), na.rm=TRUE)
    se.avg <- sd(att.avg.boot, na.rm=TRUE)
    pvalue.avg <- get.pvalue(att.avg.boot)
    est.avg <- t(as.matrix(c(att.avg, se.avg, CI.avg, pvalue.avg)))
    colnames(est.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")

    
    ## regression coefficents
    if (p > 0) {
        CI.beta<-t(apply(beta.boot, 1, function(vec)
            quantile(vec,c(0.025, 0.975), na.rm=TRUE)))
        se.beta<-apply(beta.boot, 1, function(vec)sd(vec,na.rm=TRUE))
        pvalue.beta <- apply(beta.boot, 1, get.pvalue)
        beta[na.pos] <- NA
        est.beta<-cbind(beta, se.beta, CI.beta, pvalue.beta)
        colnames(est.beta)<-c("beta", "S.E.", "CI.lower", "CI.upper", "p.value")
    }

    ## placebo test
    if (!is.null(placebo.period) & placeboTest == TRUE) {
        att.placebo <- out$att.placebo        
        CI.placebo <- quantile(att.placebo.boot, c(0.025,0.975), na.rm=TRUE)
        se.placebo <- sd(att.placebo.boot, na.rm=TRUE)
        pvalue.placebo <- get.pvalue(att.placebo.boot)
        est.placebo <- t(as.matrix(c(att.placebo, se.placebo, CI.placebo, pvalue.placebo)))
        colnames(est.placebo) <- c("ATT.placebo", "S.E.", "CI.lower", "CI.upper", "p.value")
    }
  
    ##storage
    result<-list(est.avg = est.avg,
                 att.avg.boot = att.avg.boot,
                 est.att.on = est.att.on,
                 att.on.boot = att.on.boot,
                 att.on.count.boot = att.on.count.boot,
                 T0.on.p = T0.on.p)
    if (p>0) {
        result <- c(result,list(beta.boot = beta.boot))
        result <- c(result,list(est.beta = est.beta))
    }
    if (hasRevs == 1) {
        result<-c(result,list(est.att.off = est.att.off, att.off.boot = att.off.boot, att.off.count.boot = att.off.count.boot, T0.off.p = T0.off.p))
    } 

    if (!is.null(placebo.period) & placeboTest == TRUE) {
        result <- c(result, list(est.placebo = est.placebo, att.placebo.boot = att.placebo.boot))
    }

    return(c(out,result))

    
} ## end of boot