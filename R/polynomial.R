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
                            power = 1,
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

    oci <- which(c(II) == 1)

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
    vindex <- cbind(rep(1:N, each = TT), rep(1:TT - 1, N))  ## id time
    if (power > 1) {
        for (i in 2:power) {
            vindex <- cbind(vindex, rep((1:TT - 1)^i, N))
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

    if (power > 1) {
        for (i in 2:power) {
            cf <- c(cf, list(c(1, i + 1)))
        }
    }

    est.best <- fastplm(y = as.matrix(vy[oci]), 
                        x = vx.fit, 
                        ind = as.matrix(vindex[oci,]),
                        sfe = sf, cfe = cf, PCA = TRUE,
                        se = FALSE)

    yfit <- predict(est.best, x = vx, ind = vindex)

    Y.ct <- matrix(yfit, TT, N)

    beta <- est.best$coefficients
    
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

    }

    ## 0. relevant parameters
    if (is.null(beta)) {
        beta <- NA
    }
   
    ## 1. estimated att and counterfactuals
    eff <- Y - Y.ct    
    att.avg <- sum(eff * D)/(sum(D))

    ## 2. rmse for treated units' observations under control
    tr <- which(apply(D, 2, sum) > 0)
    tr.co <- which((as.matrix(1 - D[,tr]) * as.matrix(II[,tr])) == 1)
    eff.tr <- as.matrix(eff[,tr])
    v.eff.tr <- eff.tr[tr.co]
    rmse <- sqrt(mean(v.eff.tr^2))
    

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
  
    ##-------------------------------##
    ##            Storage            ##
    ##-------------------------------##  
    out<-list(
        ## main results 
        method = "polynomial",
        Y.ct = Y.ct,
        eff = eff,
        att.avg = att.avg,
        ## supporting
        force = force,
        T = TT,
        N = N,
        p = p,
        beta = beta,
        est = est.best,
        validX = validX,
        time = time.on,
        att = att.on,
        count = count.on,
        eff.pre = eff.pre,
        eff.pre.equiv = eff.pre.equiv,
        sigma2.pre = sigma2.pre)
    
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

    return(out)
} 




