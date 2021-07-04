###################################################################
## IFE Model Function
###################################################################
fect.polynomial <- function(Y, # Outcome variable, (T*N) matrix
                            X, # Explanatory variables:  (T*N*p) array
                            D, #  Indicator for treated unit (tr==1) 
                            I,
                            II, 
                            T.on, 
                            T.off = NULL, 
                            method,
                            degree = 1,
                            sfe = NULL,
                            cfe = NULL,
                            ind.matrix = NULL,
                            knots = NULL,
                            force, 
                            hasRevs = 1,
                            tol, # tolerance level
                            boot = FALSE, # bootstrapped sample
                            placeboTest = 0,
                            placebo.period = NULL,
                            carryoverTest = 0,
                            carryover.period = NULL,
                            norm.para = NULL,
                            time.on.seq = NULL,
                            time.off.seq = NULL,
                            group.level = NULL,
                            group = NULL,
                            time.on.seq.group = NULL,
                            time.off.seq.group = NULL) {  
    
    ##-------------------------------##
    ## Parsing data
    ##-------------------------------##  
    carryover.pos <- placebo.pos <- na.pos <- NULL
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
    est.fect <- NULL
    if (boot == FALSE) {
        est.fect <- inter_fe_ub(YY, Y0, X, II, beta0, 0, force = force, tol)
    }

    ## reshape 
    vy <- as.matrix(c(YY))
    vx.fit <- vx <- NULL
    if (p > 0) {
        vx <- matrix(NA, N*TT, p)
        for (i in 1:p) {
            vx[, i] <- c(X[,, i])
        }
        vx.fit <- as.matrix(vx[oci,])
    }
    
    vindex <- NULL
    sp <- NULL
    sf <- NULL
    cf <- NULL

    if(method == "cfe"){
        vindex <- cbind(rep(1:N, each = TT), rep(1:TT, N))  ## id time
        if (force == 1) {
            sf <- 1
        } 
        else if (force == 2) {
            sf <- 2
        } 
        else {
            sf <- c(1,2)
        }

        ## simple fixed effects
        for(ind.name in names(ind.matrix)){
            vindex <- cbind(vindex,matrix(ind.matrix[[ind.name]],ncol=1))
        }

        ind.name <- c("id","time",names(ind.matrix))
        ind.index <- c(1:(2+length(names(ind.matrix))))
        names(ind.index) <- ind.name

        if(!is.null(sfe)){
            for(sub.sfe in sfe){
                sf.add <- ind.index[sub.sfe]
                names(sf.add) <- NULL
                sf <- c(sf,sf.add)
            }
        }

        if(!is.null(cfe)){
            for(sub.cfe in cfe){
                sub.cf <- c(ind.index[sub.cfe[1]], ind.index[sub.cfe[2]])
                names(sub.cf) <- NULL
                cf <- c(cf, list(sub.cf))
            }
        }

        est.best <- suppressWarnings(fastplm(y = as.matrix(vy[oci]), 
                                             x = vx.fit, 
                                             ind = as.matrix(vindex[oci,]),
                                             sfe = sf, cfe = cf, PCA = TRUE,
                                             se = FALSE,
                                             drop.singletons = FALSE))
        yfit <- predict(est.best, x = vx, ind = vindex)

    }
    else if (method == "polynomial") {
        vindex <- cbind(rep(1:N, each = TT), rep(1:TT, N))  ## id time
        if (degree > 1) {
            for (i in 2:degree) {
                vindex <- cbind(vindex, rep((1:TT)^i, N))
            }
        }
        if (force == 1) {
            sf <- 1
        } else if (force == 2) {
            sf <- 2
        } else {
            sf <- c(1,2)
        }

        cf <- list(c(1,2))

        if (degree > 1) {
            for (i in 2:degree) {
                cf <- c(cf, list(c(1, i + 1)))
            }
        }

        est.best <- suppressWarnings(fastplm(y = as.matrix(vy[oci]), 
                                            x = vx.fit, 
                                            ind = as.matrix(vindex[oci,]),
                                            sfe = sf, cfe = cf, PCA = TRUE,
                                            se = FALSE))

        yfit <- predict(est.best, x = vx, ind = vindex)


    } 
    else {
        sf <- 1
        vindex <- as.matrix(rep(1:N, each = TT))
        sp <- as.matrix(rep(1:TT, N))
        est.best <- suppressWarnings(fastplm(y = as.matrix(vy[oci]), 
                                            x = vx.fit, 
                                            ind = as.matrix(vindex[oci,]),
                                            sp = as.matrix(sp[oci,]),
                                            degree = degree,
                                            sfe = sf, cfe = cf, PCA = 0,
                                            se = FALSE))

        yfit <- predict(est.best, x = vx, ind = vindex, sp = sp)

    }


    Y.ct <- matrix(yfit, TT, N)
    if (p > 0) {
        beta <- as.matrix(c(est.best$coefficients)[1:p])
    } else {
        beta <- matrix(0, 1, 0)
    }
    est.best$beta <- beta
    validX <- ifelse(p > 0, 1, 0) 

    ##-------------------------------##
    ##   ATT and Counterfactuals     ##
    ##-------------------------------##

    ## we first adjustment for normalization 
    if (!is.null(norm.para)) {
        Y <- Y * norm.para[1]
        ## variance of the error term 
        sigma2 <- est.best$sigma2 * (norm.para[1]^2)    
        est.best$sigma2 <- sigma2
        ## output of estimates
        est.best$mu <- est.best$mu * norm.para[1] 
        est.best$residuals <- est.best$residuals * norm.para[1]
        Y.ct <- Y.ct * norm.para[1]
        if (boot == FALSE) {
            est.fect$fit <- est.fect$fit * norm.para[1]
        }
        est.fect$sigma2 <- est.fect$sigma2 * norm.para[1]
    }

    ## 0. relevant parameters
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
    eff <- Y - Y.ct   
    complete.index <- which(!is.na(eff))
    missing.index <- which(is.na(eff))
    if(length(missing.index)>0){
        I[missing.index] <- 0
        II[missing.index] <- 0
    }
    att.avg <- sum(eff[complete.index] * D[complete.index])/(sum(D[complete.index]))

    ## att.avg.unit
    tr.pos <- which(apply(D, 2, sum) > 0)
    att.unit <- sapply(1:length(tr.pos), function(vec){return((eff[, tr.pos[vec]] * D[, tr.pos[vec]]) / sum(D[, tr.pos[vec]]))})
    att.avg.unit <- mean(att.unit,na.rm=TRUE)


    ## 2. rmse for treated units' observations under control
    tr <- which(apply(D, 2, sum) > 0)
    tr.co <- which((as.matrix(1 - D[,tr]) * as.matrix(II[,tr])) == 1)
    eff.tr <- as.matrix(eff[,tr])
    v.eff.tr <- eff.tr[tr.co]
    rmse <- sqrt(mean(v.eff.tr^2,na.rm=TRUE))
    

    ## 3. unbalanced output
    if (0 %in% I) {
        eff[which(I == 0)] <- NA
        Y.ct[which(I == 0)] <- NA
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

    sigma2.pre <- eff.pre.equiv <- NULL
    if (boot == FALSE) {
        eff.pre.equiv <- eff.pre

        sigma2.pre <- tapply(eff.pre.equiv[,1], eff.pre.equiv[,2], var)
        sigma2.pre <- cbind(sigma2.pre, sort(unique(eff.pre.equiv[, 2])), table(eff.pre.equiv[, 2]))
        colnames(sigma2.pre) <- c("sigma2", "period", "count")
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

        if (boot == FALSE) {
            eff.off.equiv <- eff.off

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

    ## 7. carryover effects
    if (!is.null(carryover.period) && carryoverTest == 1 && hasRevs) {              
        if (length(carryover.period) == 1) {
            carryover.pos <- which(time.off == carryover.period)
            att.carryover <- att.off[carryover.pos]
        } else {
            carryover.pos <- which(time.off >= carryover.period[1] & time.off <= carryover.period[2])
            att.carryover <- sum(att.off[carryover.pos] * count.off[carryover.pos]) / sum(count.off[carryover.pos])
        }
    }
    
    ## 8. cohort effects
    if (!is.null(group)) {
        cohort <- cbind(c(group), c(D), c(eff.v))
        rm.pos <- unique(c(rm.pos1, which(cohort[, 2] == 0)))
        cohort <- cohort[-rm.pos, ]

        g.level <- sort(unique(cohort[, 1]))
        raw.group.att <- as.numeric(tapply(cohort[, 3], cohort[, 1], mean))

        group.att <- rep(NA, length(group.level))
        group.att[which(group.level %in% g.level)] <- raw.group.att

        # by-group dynamic effects
        group.level.name <- names(group.level)

        group.output <- list()
        for(i in c(1:length(group.level))){
            sub.group <- group.level[i]
            sub.group.name <- group.level.name[i]

            ## by-group dynamic effects
            t.on.sub <- c(T.on[which(group==sub.group)])
            eff.v.sub <- c(eff[which(group==sub.group)]) ## a vector
            rm.pos1.sub <- which(is.na(eff.v.sub))
            rm.pos2.sub <- which(is.na(t.on.sub)) 
            eff.v.use1.sub <- eff.v.sub
            t.on.use.sub <- t.on.sub
            if (NA %in% eff.v.sub | NA %in% t.on.sub) {
                eff.v.use1.sub <- eff.v.sub[-c(rm.pos1.sub, rm.pos2.sub)]
                t.on.use.sub <- t.on.sub[-c(rm.pos1.sub, rm.pos2.sub)]
            }
            if(length(t.on.use.sub)>0){
                time.on.sub <- sort(unique(t.on.use.sub))
                att.on.sub <- as.numeric(tapply(eff.v.use1.sub, 
                                            t.on.use.sub, 
                                            mean)) ## NA already removed
                count.on.sub <- as.numeric(table(t.on.use.sub))
            }else{
                time.on.sub <- att.on.sub <- count.on.sub <- NULL
            }
            
            if (!is.null(time.on.seq.group)) {
                count.on.med.sub <- att.on.med.sub <- rep(NA, length(time.on.seq.group[[sub.group.name]]))
                time.on.seq.sub <- time.on.seq.group[[sub.group.name]]
                att.on.med.sub[which(time.on.seq.sub %in% time.on.sub)] <- att.on.sub
                count.on.med.sub[which(time.on.seq.sub %in% time.on.sub)] <- count.on.sub
                att.on.sub <- att.on.med.sub
                count.on.sub <- count.on.med.sub
                time.on.sub<- time.on.seq.sub
            }
            suboutput <- list(att.on=att.on.sub,
                              time.on=time.on.sub,
                              count.on=count.on.sub)

            ## placebo effect, if placeboTest == 1 
            if (!is.null(placebo.period) && placeboTest == 1) {              
                if (length(placebo.period) == 1) {
                    placebo.pos.sub <- which(time.on.sub == placebo.period)
                    if(length(placebo.pos.sub)>0){
                        att.placebo.sub <- att.on.sub[placebo.pos.sub]
                    }
                    else{att.placebo.sub <- NULL} 
                } 
                else {
                    placebo.pos.sub <- which(time.on.sub >= placebo.period[1] & time.on.sub <= placebo.period[2])
                    if(length(placebo.pos.sub)>0){
                        att.placebo.sub <- sum(att.on.sub[placebo.pos.sub] * count.on.sub[placebo.pos.sub]) / sum(count.on.sub[placebo.pos.sub])
                    }
                    else{att.placebo.sub <- NULL} 
                }
                suboutput <- c(suboutput, list(att.placebo = att.placebo.sub))
            }

            ## T.off
            if (hasRevs == 1) {    
                t.off.sub <- c(T.off[which(group==sub.group)])
                rm.pos3.sub <- which(is.na(t.off.sub))
                eff.v.use2.sub <- eff.v.sub
                t.off.use.sub <- t.off.sub
                if (NA %in% eff.v.sub | NA %in% t.off.sub) {
                    eff.v.use2.sub <- eff.v.sub[-c(rm.pos1.sub, rm.pos3.sub)]
                    t.off.use.sub <- t.off.sub[-c(rm.pos1.sub, rm.pos3.sub)]
                }
                if(length(t.off.use.sub)>0){
                    time.off.sub <- sort(unique(t.off.use.sub))
                    att.off.sub <- as.numeric(tapply(eff.v.use2.sub, t.off.use.sub, mean)) ## NA already removed
                    count.off.sub <- as.numeric(table(t.off.use.sub))
                }else{
                    time.off.sub <- att.off.sub <- count.off.sub <- NULL
                }

                if (!is.null(time.off.seq.group)) {
                    count.off.med.sub <- att.off.med.sub <- rep(NA, length(time.off.seq.group[[sub.group.name]]))
                    time.off.seq.sub <- time.off.seq.group[[sub.group.name]]
                    att.off.med.sub[which(time.off.seq.sub %in% time.off.sub)] <- att.off.sub
                    count.off.med.sub[which(time.off.seq.sub %in% time.off.sub)] <- count.off.sub
                    att.off.sub <- att.off.med.sub
                    count.off.sub <- count.off.med.sub
                    time.off.sub <- time.off.seq.sub
                }
                suboutput <- c(suboutput, list(att.off = att.off.sub,
                                               count.off = count.off.sub,
                                               time.off = time.off.sub))

                if (!is.null(carryover.period) && carryoverTest == 1) {
                    if (length(carryover.period) == 1) {
                        carryover.pos.sub <- which(time.off.sub == carryover.period.sub)
                        if(length(carryover.pos.sub)>0){
                            att.carryover.sub <- att.off.sub[carryover.pos.sub]
                        } else{att.carryover.sub <- NULL}      
                    } else {
                        carryover.pos.sub <- which(time.off.sub >= carryover.period[1] & time.off.sub <= carryover.period[2])
                        if(length(carryover.pos.sub)>0){
                            att.carryover.sub <- sum(att.off.sub[carryover.pos.sub] * count.off.sub[carryover.pos.sub]) / sum(count.off.sub[carryover.pos.sub])
                        } else{att.carryover.sub <- NULL}  
                    }
                    suboutput <- c(suboutput,list(att.carryover = att.carryover.sub))
                }
            }
            group.output[[sub.group.name]] <- suboutput
        }
    }

  
    ##-------------------------------##
    ##            Storage            ##
    ##-------------------------------##  
    out<-list(
        ## main results 
        method = method,
        Y.ct = Y.ct,
        eff = eff,
        I = I,
        II = II,
        att.avg = att.avg,
        att.avg.unit = att.avg.unit,
        ## supporting
        force = force,
        T = TT,
        N = N,
        p = p,
        beta = beta,
        est = est.best,
        sigma2 = est.best$sigma2,
        sigma2.fect = est.fect$sigma2,
        validX = validX,
        time = time.on,
        att = att.on,
        count = count.on,
        eff.pre = eff.pre,
        eff.pre.equiv = eff.pre.equiv,
        sigma2.pre = sigma2.pre)
    
    #print(att.on)

    if (hasRevs == 1) {
        out <- c(out, list(time.off = time.off, 
                           att.off = att.off,
                           count.off = count.off,
                           eff.off = eff.off,
                           eff.off.equiv = eff.off.equiv,
                           off.sd = off.sd))
    }

    if (!is.null(placebo.period) && placeboTest == 1) {
        out <- c(out, list(att.placebo = att.placebo))
    }

    if (!is.null(carryover.period) && carryoverTest == 1) {
        out <- c(out, list(att.carryover = att.carryover))
    }

    if (!is.null(group)) {
        out <- c(out, list(group.att = group.att,
                           group.output = group.output))
    }
    return(out)
} 




