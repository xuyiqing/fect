#################################
## support function
#################################
get_term <- function(d,
                     ii, 
                     type = "on") {
    dd <- d
    iii <- ii
    ## dd <- dd[which(iii == 1)]
    first.pos <- min(which(iii == 1))
    if (first.pos != 1) {
        dd <- dd[-(1:(first.pos-1))]
        iii <- iii[-(1:(first.pos-1))]
    } 
    T <- length(dd) 
    if (0 %in% iii) {
        if (T > 1) {
            for (i in 1:(T-1)) {
                if (iii[i+1] == 0) {
                    dd[i+1] <- dd[i]
                }
            }
        }
    }

    if (type == "off") {
        dd <- abs(dd - 1)
    }
    d1 <- dd[1:(T-1)]
    d2 <- dd[2:T]
    
    if(T==1){
        term <- rep(NA, 1)
    }
    else if (sum(d1 == d2) == (T-1)) {
        term <- rep(NA, T)
    } 
    else {
        change.pos <- which(d1 != d2) + 1
        change.length <- length(change.pos)
        term <- NULL    
        if (dd[1] == 0) {   
            for (i in 1:(change.length)) {
                if (i == 1) {
                    part.term <- (2 - change.pos[i]):0
                } else {
                    if (i %% 2 == 0) {
                        part.term <- 1:(change.pos[i] - change.pos[i-1])
                    } else {
                        part.term <- (change.pos[i-1] - change.pos[i] + 1):0
                    }
                }
                term <- c(term, part.term)
            }
        } else if (dd[1] == 1) {
            for (i in 1:(change.length)) {
                if (i == 1) {
                    #if (type == "on") {
                    #    part.term <- 1:(change.pos[i] - 1)
                    #} else if (type == "off") {
                        part.term <- rep(NA, change.pos[i] - 1)
                    #}
                } else {
                    if (i %% 2 == 0) {
                        part.term <- (change.pos[i-1] - change.pos[i] + 1):0
                    } else {
                        part.term <- 1:(change.pos[i] - change.pos[i-1])
                    }
                }
                term <- c(term, part.term)
            }
        } 
        if (dd[change.pos[change.length]] == 0) {
            term <- c(term, rep(NA, (T - change.pos[change.length] + 1)))
        } else {
            term <- c(term, 1:(T - change.pos[change.length] + 1))
        }
    }
    ## term.all <- rep(NA, length(d))
    if (first.pos != 1) {
        term <- c(rep(NA, (first.pos - 1)), term)
    }
    return(term)
}

###################################
## regressions for initial values
###################################

initialFit <- function(data, ## long form data 
                       force,
                       oci ) { ## indicator

    N <- length(unique(data[,2]))
    T <- length(unique(data[,3]))
    p <- dim(data)[2] - 3

    x <- x.sub <- NULL
    if (p > 0) {
        x <- as.matrix(data[,4:(dim(data)[2])])
        x.sub <- as.matrix(x[oci,])
    }
    y <- as.matrix(data[,1])
    beta0 <- matrix(0, 1, 1)

    ind <- NULL
    if(force == 0){
        ind <- NULL
    }
    else if (force == 1) {
        ind <- as.matrix(data[,2])
    } else if (force == 2) {
        ind <- as.matrix(data[,3])
    } else if (force == 3) {
        ind <- as.matrix(data[,c(2,3)])
    }

    if (force == 0) {
        if (p == 0) {
            mu <- mean(c(y)[oci])
            Y0 <- matrix(mu, T, N)
            ## res <- as.matrix(c(y) - mu)
        } 
        else {
            lm.fit <- lm(as.matrix(c(y)[oci])~x.sub)
            coef <- lm.fit$coefficients
            mu <- coef[1]
            beta0 <- as.matrix(coef[2:length(coef)])
            if (sum(is.na(beta0)) > 0) {
                beta0[which(is.na(beta0))] <- 0
            }
            y0 <- mu + x %*% beta0
            Y0 <- matrix(y0, T, N)
            ## res <- as.matrix(lm.fit$residuals)
        }
    } 
    else {
        
        colnames(y) <- y.name <- 'y'
        colnames(ind) <- ind.name <- paste0("id.",c(1:dim(ind)[2]))
        if(p>0){
            colnames(x) <- x.name <- paste0("x.",c(1:dim(x)[2]))
            data.reg <- cbind.data.frame(y,x,ind)
            formula.reg <- paste0("y~",paste(x.name,collapse="+"),"|",paste(ind.name,collapse="+"))
        }
        else{
            data.reg <- cbind.data.frame(y,ind)
            formula.reg <- paste0("y~1|",paste(ind.name,collapse="+"))
        }
        formula.reg <- as.formula(formula.reg)


        lm.fit <- suppressWarnings(invisible(feols(fml = formula.reg,
                                                   data = data.reg[oci,],
                                                   fixef.rm = "none")))
        
        y0 <- suppressWarnings(predict(lm.fit, newdata = data.reg))
        #lm.fit <- suppressWarnings(invisible(fastplm(y = as.matrix(c(y)[oci]), x = x.sub, 
        #                   ind = as.matrix(ind[oci,]),drop.singletons = FALSE)))
        #y0 <- suppressWarnings(predict(lm.fit, x = x, ind = ind))
        Y0 <- matrix(y0, T, N)
        if (p > 0) {
            beta0 <- lm.fit$coefficients
            names(beta0) <- NULL
            beta0 <- as.matrix(beta0)
        }
    }
    result <- list(Y0 = Y0, beta0 = beta0)
    return(result)
}

################################################
##  regressions for initial values, probit model  ##
################################################

## if we do QR : the fixed effects term doesn't contain mu, and we need Y0, FE0, xi0, factor0
## if we do SVD : the fixed effects term contains mu, and we need Y0, FE0

BiInitialFit <- function(data, ## long form data 
                         QR = 0,
                         force,
                         r,
                         oci ) { ## indicator

    N <- length(unique(data[,2]))
    T <- length(unique(data[,3]))
    p <- dim(data)[2] - 3

    x <- x.sub <- NULL
    if (p > 0) {
        x <- as.matrix(data[,4:(dim(data)[2])])
        ## regard as missing 
        x[setdiff(1:(N*T), oci),] <- 0
        x.sub <- as.matrix(x[oci,])
    }
    y <- as.matrix(data[,1])
    beta0 <- matrix(0, 1, 1)

    ind <- NULL
    if (force == 1) {
        ind <- as.matrix(data[,2])
    } else if (force == 2) {
        ind <- as.matrix(data[,3])
    } else if (force == 3) {
        ind <- as.matrix(data[,c(2,3)])
    }

    xi <- matrix(0, T, 1)
    alpha <- matrix(0, N, 1) 
    factor <- matrix(0, T, r)
    loadings <- matrix(0, N, r)
    Y0 <- FE <- matrix(0, T, N)

    if (force == 0) { ## no additive fixed effects
        if (p == 0) {
            mu <- mean(c(y)[oci])
            y0 <- as.matrix(rep(mu, T*N))
            Y0 <- matrix(mu, T, N)
        } else {
            lm.fit <- lm(as.matrix(c(y)[oci])~x.sub)
            coef <- lm.fit$coefficients
            mu <- coef[1]
            beta0 <- as.matrix(coef[2:length(coef)])
            if (sum(is.na(beta0)) > 0) {
                beta0[which(is.na(beta0))] <- 0
            }
            y0 <- mu + x %*% beta0
            Y0 <- matrix(y0, T, N)
        }
        if (QR == 0) {
            FE <- mu
        }
    } else {         ## with additive fixed effects
        #plm.fit <- suppressWarnings(invisible(fastplm(y = as.matrix(c(y)[oci]), x = x.sub, 
        #                   ind = as.matrix(ind[oci,]),drop.singletons = FALSE)))
        plm.fit <- NULL #to delete
        y0 <- suppressWarnings(predict(plm.fit, x = x, ind = ind))
        Y0 <- matrix(y0, T, N)
        mu <- plm.fit$intercept
        if (p > 0) {
            beta0 <- plm.fit$coefficients
        }
        if (force == 1) {
            alpha <- plm.fit$sfe.coefs[[1]]
        }
        if (force == 2) {
            xi <- plm.fit$sfe.coefs[[1]]
        } 
        else if (force == 3) {
            alpha <- plm.fit$sfe.coefs[[1]]
            xi <- plm.fit$sfe.coefs[[2]]
        } 
    }

    ## pca 
    if (r > 0) { ## factor analysis of residuals
        
        res <- y - y0
        res[setdiff(1:(N*T), oci),] <- 0 
        res <- matrix(res, T, N)        
        ife_pca <- panel_factor(res, r) 
        ife <- ife_pca$FE 
        factor <- ife_pca$factor
        loadings <- ife_pca$lambda
        Y0 <- Y0 + ife 

        if (QR == 1) {
            ## initial value: need time fixed effects xi and factor f
            if (force == 0) {
                ## factor <- qr_factor(factor, loadings)$factor
                FE <- ife_pca$FE 
            }
            else if (force == 1) { 
                ## modify initial Y0 
                Y0[setdiff(1:(N*T), oci)] <- Y0[setdiff(1:(N*T), oci)] - (sum(factor[1,] * loadings[1,]) + alpha[1])
                ## restrctions: alpha_1 = 0, f_1 = 0 
                alpha <- alpha - alpha[1] + (loadings - matrix(rep(loadings[1, ], N), N, r, byrow = TRUE)) %*% as.matrix(factor[1, ])
                factor <- factor - matrix(rep(factor[1, ], T), T, r, byrow = TRUE)
                ife_qr <- qr_factor(factor, loadings)
                factor <- ife_qr$factor 
                FE <- ife_qr$FE + matrix(rep(alpha, each = T), T, N)
            }
            else if (force == 2) { 
                ## modify initial Y0 
                Y0[setdiff(1:(N*T), oci)] <- Y0[setdiff(1:(N*T), oci)] - (sum(factor[1,] * loadings[1,]) + xi[1])
                ## restrctions: xi_1 = 0 
                xi <- xi - xi[1] + (factor - matrix(rep(factor[1, ], T), T, r, byrow = TRUE)) %*% as.matrix(loadings[1, ])
                loadings <- loadings - matrix(rep(loadings[1, ], N), N, r, byrow = TRUE)
                ife_qr <- qr_factor(factor, loadings)
                factor <- ife_qr$factor
                FE <- ife_qr$FE + matrix(rep(xi, N), T, N)
            }
            else if (force == 3) {
                ## modify initial Y0 
                Y0[setdiff(1:(N*T), oci)] <- Y0[setdiff(1:(N*T), oci)] - (sum(factor[1,] * loadings[1,]) + xi[1] + alpha[1])
                ## restrctions: xi_1 = 0 , f_1 = 0 
                alpha <- alpha - alpha[1] + (loadings - matrix(rep(loadings[1, ], N), N, r, byrow = TRUE)) %*% as.matrix(factor[1, ])
                xi <- xi - xi[1] + (factor - matrix(rep(factor[1, ], T), T, r, byrow = TRUE)) %*% as.matrix(loadings[1, ])
                factor <- factor - matrix(rep(factor[1, ], T), T, r, byrow = TRUE)
                loadings <- loadings - matrix(rep(loadings[1, ], N), N, r, byrow = TRUE)
                ife_qr <- qr_factor(factor, loadings)
                factor <- ife_qr$factor
                FE <- ife_qr$FE + matrix(rep(alpha, each = T), T, N) + matrix(rep(xi, N), T, N)
            }
        } else {
            FE <- mu + matrix(rep(alpha, each = T), T, N) + matrix(rep(xi, N), T, N) + ife
        }

    } else { ## only adjust additive fixed effects 
        if (QR == 1) {
            if (force == 1) {
                ## modify initial Y0 
                Y0[setdiff(1:(N*T), oci)] <- Y0[setdiff(1:(N*T), oci)] - alpha[1]
                ## restrctions: alpha_1 = 0 
                FE <- matrix(rep(alpha, each = T), T, N)
            }
            else if (force == 2) {
                ## modify initial Y0 
                Y0[setdiff(1:(N*T), oci)] <- Y0[setdiff(1:(N*T), oci)] - xi[1]
                ## restrctions: xi_1 = 0 
                xi <- xi - xi[1] 
                FE <- matrix(rep(xi, N), T, N)
            }
            else if (force == 3) {
                ## modify initial Y0 
                Y0[setdiff(1:(N*T), oci)] <- Y0[setdiff(1:(N*T), oci)] - (xi[1] + alpha[1]) 
                ## restrctions: xi_1 = 0 , alpha_1 = 0
                xi <- xi - xi[1] 
                FE <- matrix(rep(alpha, each = T), T, N) + matrix(rep(xi, N), T, N)
            }
        } else {
            FE <- mu + matrix(rep(alpha, each = T), T, N) + matrix(rep(xi, N), T, N)
        }
        
    }

    result <- list(Y0 = Y0, FE0 = FE, xi0 = xi, factor0 = factor, loadings = loadings)
    return(result)
}

## cross validation sampling
cv.sample2 <- function(I, count) {
    N <- dim(I)[2]
    TT <- dim(I)[1]
    cv.id <- NULL
    oci <- which(c(I) == 1)
    if (count <= 3) {
        cv.id <- sample(oci, count, replace = FALSE)   
    } else {
        ## remove boundary observation 
        oci2 <- setdiff(oci, c((1:N)*TT, (TT*(0:(N-1))+1)))
        ## randomly select 1/3
        subcount <- floor(count/3)
        rm.id <- sample(oci2, subcount, replace = FALSE)
        rm.id.upper <- rm.id - 1
        rm.id.lower <- rm.id + 1

        cv.id1 <- unique(c(rm.id, rm.id.upper, rm.id.lower))
        pos <- sapply(1:length(cv.id1), function(i) cv.id1[i]%in%oci)
        cv.id1 <- cv.id1[pos]
        if (length(cv.id1) >= count) {
            cv.id1 <- cv.id1[1:count]
            cv.id2 <- NULL
        } else {
            cv.id2 <- sample(setdiff(oci, cv.id1), (count - length(cv.id1)), replace = FALSE)
        }
        cv.id <- sort(c(cv.id1, cv.id2))
    }
    return(cv.id)
}



## cross validation sampling
#cv.sample <- function(I, D, count, 
#                      cv.count = 3, 
#                      cv.treat = FALSE) {
    
#    N <- dim(I)[2]
#    TT <- dim(I)[1]
#    tr.pos <- which(apply(D, 2, sum) >= 1) ## treated units
#    D.fake <- matrix(0, TT, N)

#    cv.id <- NULL
#    if (cv.treat == FALSE) {
#        oci <- which(c(I) == 1)
#    } else {
#        D.fake[, tr.pos] <- 1
#        oci <- which(c(I) == 1 & c(D.fake) == 1)
#    }
    
#    if (length(oci) <= count) {
#        stop("Too few observations are valid for cross-validation. Try setting the option cv.treat to FALSE.\n")
#    }

#    if (cv.count == 1 || count <= 2) {  ## randomly missing
#        cv.id <- sample(oci, count, replace = FALSE)
#    } else {
        ## remove boundary observation 
#        if (cv.treat == FALSE) {
#            rm.pos <- c()
#            for (i in 1:(cv.count - 1)) {
#                rm.pos <- c(rm.pos, (TT * (0:(N-1)) + i))
#            }
#            oci2 <- setdiff(oci, rm.pos)
#        } else {
#            rm.pos <- c()
#            for (i in 1:(cv.count - 1)) {
#                rm.pos <- c(rm.pos, (TT * (tr.pos - 1) + i))
#            }
#            oci2 <- setdiff(oci, rm.pos)
#        }

        ## randomly select 1/cv.count
#        subcount <- floor(count/cv.count)
#        if (subcount == 0) {
#            subcount <- 1
#        }
#        rm.id <- sample(oci2, subcount, replace = FALSE)
#        rm.id.all <- c()
#        for (i in 1:(cv.count - 1)) {
#            rm.id.all <- c(rm.id.all, rm.id - i)
#        }
#        rm.id.all <- c(rm.id.all, rm.id)

#        cv.id1 <- unique(rm.id.all)
#        pos <- unlist(sapply(1:length(cv.id1), function(i) cv.id1[i] %in% oci))
#        cv.id1 <- cv.id1[pos]

#        if (length(cv.id1) >= count) {
#            cv.id1 <- cv.id1[1:count]
#            cv.id2 <- NULL
#        } else {
#            cv.id2 <- sample(setdiff(oci, cv.id1), (count - length(cv.id1)), replace = FALSE)
#        }
#        
#        cv.id <- sort(c(cv.id1, cv.id2))
#    }
    
#    return(cv.id)
#}

## cross validation sampling
cv.sample <- function(I, D, count, 
                      cv.count = 3, 
                      cv.treat = FALSE,
                      cv.donut = 1) {
    
    ## prop <- sum(c(I))
    N <- dim(I)[2]
    TT <- dim(I)[1]
    tr.pos <- which(apply(D, 2, sum) >= 1) ## treated units
    D.fake <- matrix(0, TT, N)

    prop <- (TT * N) / sum(c(I))

    cv.id <- NULL
    if (cv.treat == FALSE) {
        oci <- which(c(I) == 1)
    } else {
        D.fake[, tr.pos] <- 1
        oci <- which(c(I) == 1 & c(D.fake) == 1)
        if (length(oci) <= count) {
            stop("Too few observations are valid for cross-validation. Try to set the option cv.treat to FALSE or set a smaller cv.prop.\n")
        }
    }

    if (cv.count == 1 || count <= 2) {  ## randomly missing
        cv.id <- sample(oci, count, replace = FALSE)
        rm.id.use <- cv.id
    } 
    else {
        res <- TT %% cv.count 
        int <- floor(TT / cv.count)
        rm.pos <- c()

        ## randomly select 1/cv.count
        subcount <- floor(count/cv.count * prop)
        
        if (subcount == 0) {
            subcount <- 1
        }

        if (cv.treat == FALSE) {
            for (i in 1:N) {
                rm.pos <- c(rm.pos, (TT * (i - 1) + res) + seq(from = 1, by = cv.count, length.out = int))
            }
        } else {
            for (i in tr.pos) {
                rm.pos <- c(rm.pos, (TT * (i - 1) + res) + seq(from = 1, by = cv.count, length.out = int))
            }
        }

        rm.id <- sample(rm.pos, subcount, replace = FALSE)
        rm.id.all <- rm.id 
        rm.id.use <- NULL

        if (cv.count == 2) {
            rm.id.all <- c(rm.id, rm.id + 1)
            rm.id.use <- rm.id.all
        } else {
            for (i in 0:(cv.count -  1)) {
                rm.id.all <- c(rm.id.all, rm.id + i)
                if (i >= cv.donut && i <= (cv.count - cv.donut - 1)) {
                    rm.id.use <- c(rm.id.use, rm.id + i)
                }
            }
        }

        ##rm.id.all <- c(rm.id, rm.id + 1, rm.id + 2)
        
        rm.id.all <- intersect(rm.id.all, oci) ## remove missing values
        rm.id.all <- sort(rm.id.all)

        rm.id.use <- intersect(rm.id.use, oci)
        rm.id.use <- sort(rm.id.use)

        cv.id2 <- NULL
        if (length(rm.id.all) >= count) {
            cv.id <- rm.id.all[1:count]
            pos.cv <- which(rm.id.use >= min(cv.id) & rm.id.use <= max(cv.id))
            rm.id.use <- rm.id.use[pos.cv]
            #if (length(rm.id.use) >= count) {
            #    rm.id.use <- rm.id.use[1:count]
            #}
        } else {
            cv.id2 <- sample(setdiff(oci, rm.id.all), (count - length(rm.id.all)), replace = FALSE)
            cv.id <- sort(c(rm.id.all, cv.id2))
            ## rm.id.use <- sort(c(rm.id.use, cv.id2))

        }

        ## remove boundary observation 
        #if (cv.treat == FALSE) {
        #    rm.pos <- c()
        #    for (i in 1:(cv.count - 1)) {
        #        rm.pos <- c(rm.pos, (TT * (0:(N-1)) + i))
        #    }
        #    oci2 <- setdiff(oci, rm.pos)
        #} else {
        #    rm.pos <- c()
        #    for (i in 1:(cv.count - 1)) {
        #        rm.pos <- c(rm.pos, (TT * (tr.pos - 1) + i))
        #    }
        #    oci2 <- setdiff(oci, rm.pos)
        #}

        ## randomly select 1/cv.count
        #subcount <- floor(count/cv.count)
        #if (subcount == 0) {
        #    subcount <- 1
        #}
        #rm.id <- sample(oci2, subcount, replace = FALSE)
        #rm.id.all <- c()
        #for (i in 1:(cv.count - 1)) {
        #    rm.id.all <- c(rm.id.all, rm.id - i)
        #}
        #rm.id.all <- c(rm.id.all, rm.id)

        #cv.id1 <- unique(rm.id.all)
        #pos <- unlist(sapply(1:length(cv.id1), function(i) cv.id1[i] %in% oci))
        #cv.id1 <- cv.id1[pos]

        #if (length(cv.id1) >= count) {
        #    cv.id1 <- cv.id1[1:count]
        #    cv.id2 <- NULL
        #} else {
        #    cv.id2 <- sample(setdiff(oci, cv.id1), (count - length(cv.id1)), replace = FALSE)
        #}
        
        #cv.id <- sort(c(cv.id1, cv.id2))
    }
    
    return(list(cv.id = cv.id,         ## marked id
                est.id = rm.id.use))   ## id used for mspe
}

res.vcov <- function(res, ## TT*Nboots
                     cov.ar = 1) {
    T <- dim(res)[1]
    I <- is.na(res)
    count <- matrix(NA,T,T)

    res[is.na(res)] <- 0
    vcov <- res%*%t(res)


    for (i in 1:T) {
        for (j in 1:T) {
            if (i > j) {
                count[i, j] <- count[j, i]
            } else {
                if ((j-i) <= cov.ar) {
                  II <- I[i,] + I[j,]
                  count[i, j] <- min( 1/sum(II==0), 1)
                } else {
                  count[i, j] <- 0
                }
            }
        }
    }
    vcov <- vcov*count
    return(vcov)
}


v_replace <- function(needle, haystack) { 
  sieved <- which(haystack == needle[1L]) 
  for(i in seq.int(1L, length(needle) - 1L)) {
    sieved <- sieved[haystack[sieved + i] == needle[i + 1L]]
  }
  out <- rep(NA,length(haystack))
  if(length(sieved)==0){
    return(out)
  }
  for(index in sieved){
    if(!is.na(index)){
      out[c(index:(index+length(needle)-1))] <- needle
    }
  }
  return(out)
}