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
        dd <- dd[-(1:(first.pos - 1))]
        iii <- iii[-(1:(first.pos - 1))]
    }
    T <- length(dd)
    if (0 %in% iii) {
        if (T > 1) {
            for (i in 1:(T - 1)) {
                if (iii[i + 1] == 0) {
                    dd[i + 1] <- dd[i]
                }
            }
        }
    }

    if (type == "off") {
        dd <- abs(dd - 1)
    }
    d1 <- dd[1:(T - 1)]
    d2 <- dd[2:T]

    if (T == 1) {
        term <- rep(NA, 1)
    } else if (sum(d1 == d2) == (T - 1)) {
        term <- rep(NA, T)
    } else {
        change.pos <- which(d1 != d2) + 1
        change.length <- length(change.pos)
        term <- NULL
        if (dd[1] == 0) {
            for (i in 1:(change.length)) {
                if (i == 1) {
                    part.term <- (2 - change.pos[i]):0
                } else {
                    if (i %% 2 == 0) {
                        part.term <- 1:(change.pos[i] - change.pos[i - 1])
                    } else {
                        part.term <- (change.pos[i - 1] - change.pos[i] + 1):0
                    }
                }
                term <- c(term, part.term)
            }
        } else if (dd[1] == 1) {
            for (i in 1:(change.length)) {
                if (i == 1) {
                    # if (type == "on") {
                    #    part.term <- 1:(change.pos[i] - 1)
                    # } else if (type == "off") {
                    part.term <- rep(NA, change.pos[i] - 1)
                    # }
                } else {
                    if (i %% 2 == 0) {
                        part.term <- (change.pos[i - 1] - change.pos[i] + 1):0
                    } else {
                        part.term <- 1:(change.pos[i] - change.pos[i - 1])
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

align_beta0 <- function(beta_raw, p) {
    if (p <= 0) {
        return(matrix(0, 1, 0))
    }

    beta_raw <- as.numeric(beta_raw)
    beta_names <- names(beta_raw)
    beta_full <- rep(0, p)
    filled <- rep(FALSE, p)
    used_raw <- rep(FALSE, length(beta_raw))

    if (!is.null(beta_names)) {
        for (j in seq_along(beta_raw)) {
            nm <- beta_names[j]
            idx <- suppressWarnings(as.integer(sub(".*?(\\d+)$", "\\1", nm)))
            if (!is.na(idx) && idx >= 1 && idx <= p && !used_raw[j]) {
                beta_full[idx] <- beta_raw[j]
                filled[idx] <- TRUE
                used_raw[j] <- TRUE
            }
        }
    }

    remaining_pos <- which(!filled)
    remaining_raw <- which(!used_raw)
    if (length(remaining_pos) > 0 && length(remaining_raw) > 0) {
        nfill <- min(length(remaining_pos), length(remaining_raw))
        beta_full[remaining_pos[seq_len(nfill)]] <- beta_raw[remaining_raw[seq_len(nfill)]]
    }

    beta_full[which(is.na(beta_full))] <- 0
    return(as.matrix(beta_full))
}

initialFit <- function(data, ## long form data
                       force,
                       w = NULL,
                       oci) { ## indicator

    N <- length(unique(data[, 2]))
    T <- length(unique(data[, 3]))
    p <- dim(data)[2] - 3

    x <- x.sub <- NULL
    if (p > 0) {
        x <- as.matrix(data[, 4:(dim(data)[2])])
        x.sub <- as.matrix(x[oci, ])
    }
    y <- as.matrix(data[, 1])
    y.sub <- as.matrix(y[oci, ])
    beta0 <- matrix(0, 1, 1)

    if (!is.null(w)) {
        use_weight <- 1
        w <- as.matrix(w)
        w.sub <- as.matrix(w[oci, ])
    } else {
        use_weight <- 0
    }

    ind <- NULL
    if (force == 0) {
        ind <- NULL
    } else if (force == 1) {
        ind <- as.matrix(data[, 2])
    } else if (force == 2) {
        ind <- as.matrix(data[, 3])
    } else if (force == 3) {
        ind <- as.matrix(data[, c(2, 3)])
    }

    if (force == 0) {
        if (p == 0) {
            if (use_weight == 1) {
                mu <- sum(y.sub * w.sub) / sum(w.sub)
                Y0 <- matrix(mu, T, N)
            } else {
                mu <- mean(c(y)[oci])
                Y0 <- matrix(mu, T, N)
            }
            ## res <- as.matrix(c(y) - mu)
        } else {
            if (use_weight == 1) {
                lm.fit <- lm(as.matrix(c(y)[oci]) ~ x.sub, weights = w.sub)
            } else {
                lm.fit <- lm(as.matrix(c(y)[oci]) ~ x.sub)
            }
            coef <- lm.fit$coefficients
            mu <- coef[1]
            beta0 <- align_beta0(coef[2:length(coef)], p)
            y0 <- mu + x %*% beta0
            Y0 <- matrix(y0, T, N)
            ## res <- as.matrix(lm.fit$residuals)
        }
    } else {
        colnames(y) <- y.name <- "y"
        colnames(ind) <- ind.name <- paste0("id.", c(1:dim(ind)[2]))
        if (p > 0) {
            colnames(x) <- x.name <- paste0("x.", c(1:dim(x)[2]))
            data.reg <- cbind.data.frame(y, x, ind)
            formula.reg <- paste0("y~", paste(x.name, collapse = "+"), "|", paste(ind.name, collapse = "+"))
        } else {
            data.reg <- cbind.data.frame(y, ind)
            formula.reg <- paste0("y~1|", paste(ind.name, collapse = "+"))
        }
        formula.reg <- as.formula(formula.reg)

        lm.fit <- feols(fml = formula.reg, data = data.reg[oci, ], fixef.rm = "none")
        # if (use_weight == 1) {
        #     lm.fit <- suppressWarnings(invisible(feols(
        #         fml = formula.reg,
        #         data = data.reg[oci, ], weights = w.sub,
        #         fixef.rm = "none"
        #     )))
        # } else {
        #     lm.fit <- suppressWarnings(invisible(feols(
        #         fml = formula.reg,
        #         data = data.reg[oci, ],
        #         fixef.rm = "none"
        #     )))
        # }

        y0 <- predict(lm.fit, newdata = data.reg)

        # lm.fit <- suppressWarnings(invisible(fastplm(y = as.matrix(c(y)[oci]), x = x.sub,
        #                   ind = as.matrix(ind[oci,]),drop.singletons = FALSE)))
        # y0 <- suppressWarnings(predict(lm.fit, x = x, ind = ind))
        Y0 <- matrix(y0, T, N)
        if (p > 0) {
            beta0 <- align_beta0(lm.fit$coefficients, p)
        }
    }

    bad.Y0 <- !is.finite(Y0)
    if (sum(bad.Y0) > 0) {
        mu.fill <- mean(as.numeric(y.sub), na.rm = TRUE)
        if (!is.finite(mu.fill)) {
            mu.fill <- 0
        }
        Y0[bad.Y0] <- mu.fill
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
                         oci) { ## indicator

    N <- length(unique(data[, 2]))
    T <- length(unique(data[, 3]))
    r <- min(r, T, N)
    p <- dim(data)[2] - 3

    x <- x.sub <- NULL
    if (p > 0) {
        x <- as.matrix(data[, 4:(dim(data)[2])])
        ## regard as missing
        x[setdiff(1:(N * T), oci), ] <- 0
        x.sub <- as.matrix(x[oci, ])
    }
    y <- as.matrix(data[, 1])
    beta0 <- matrix(0, 1, 1)

    ind <- NULL
    if (force == 1) {
        ind <- as.matrix(data[, 2])
    } else if (force == 2) {
        ind <- as.matrix(data[, 3])
    } else if (force == 3) {
        ind <- as.matrix(data[, c(2, 3)])
    }

    xi <- matrix(0, T, 1)
    alpha <- matrix(0, N, 1)
    factor <- matrix(0, T, r)
    loadings <- matrix(0, N, r)
    Y0 <- FE <- matrix(0, T, N)

    if (force == 0) { ## no additive fixed effects
        if (p == 0) {
            mu <- mean(c(y)[oci])
            y0 <- as.matrix(rep(mu, T * N))
            Y0 <- matrix(mu, T, N)
        } else {
            lm.fit <- lm(as.matrix(c(y)[oci]) ~ x.sub)
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
    } else { ## with additive fixed effects
        # plm.fit <- suppressWarnings(invisible(fastplm(y = as.matrix(c(y)[oci]), x = x.sub,
        #                   ind = as.matrix(ind[oci,]),drop.singletons = FALSE)))
        plm.fit <- NULL # to delete
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
        } else if (force == 3) {
            alpha <- plm.fit$sfe.coefs[[1]]
            xi <- plm.fit$sfe.coefs[[2]]
        }
    }

    ## pca
    if (r > 0) { ## factor analysis of residuals

        res <- y - y0
        res[setdiff(1:(N * T), oci), ] <- 0
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
            } else if (force == 1) {
                ## modify initial Y0
                Y0[setdiff(1:(N * T), oci)] <- Y0[setdiff(1:(N * T), oci)] - (sum(factor[1, ] * loadings[1, ]) + alpha[1])
                ## restrctions: alpha_1 = 0, f_1 = 0
                alpha <- alpha - alpha[1] + (loadings - matrix(rep(loadings[1, ], N), N, r, byrow = TRUE)) %*% as.matrix(factor[1, ])
                factor <- factor - matrix(rep(factor[1, ], T), T, r, byrow = TRUE)
                ife_qr <- qr_factor(factor, loadings)
                factor <- ife_qr$factor
                FE <- ife_qr$FE + matrix(rep(alpha, each = T), T, N)
            } else if (force == 2) {
                ## modify initial Y0
                Y0[setdiff(1:(N * T), oci)] <- Y0[setdiff(1:(N * T), oci)] - (sum(factor[1, ] * loadings[1, ]) + xi[1])
                ## restrctions: xi_1 = 0
                xi <- xi - xi[1] + (factor - matrix(rep(factor[1, ], T), T, r, byrow = TRUE)) %*% as.matrix(loadings[1, ])
                loadings <- loadings - matrix(rep(loadings[1, ], N), N, r, byrow = TRUE)
                ife_qr <- qr_factor(factor, loadings)
                factor <- ife_qr$factor
                FE <- ife_qr$FE + matrix(rep(xi, N), T, N)
            } else if (force == 3) {
                ## modify initial Y0
                Y0[setdiff(1:(N * T), oci)] <- Y0[setdiff(1:(N * T), oci)] - (sum(factor[1, ] * loadings[1, ]) + xi[1] + alpha[1])
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
                Y0[setdiff(1:(N * T), oci)] <- Y0[setdiff(1:(N * T), oci)] - alpha[1]
                ## restrctions: alpha_1 = 0
                FE <- matrix(rep(alpha, each = T), T, N)
            } else if (force == 2) {
                ## modify initial Y0
                Y0[setdiff(1:(N * T), oci)] <- Y0[setdiff(1:(N * T), oci)] - xi[1]
                ## restrctions: xi_1 = 0
                xi <- xi - xi[1]
                FE <- matrix(rep(xi, N), T, N)
            } else if (force == 3) {
                ## modify initial Y0
                Y0[setdiff(1:(N * T), oci)] <- Y0[setdiff(1:(N * T), oci)] - (xi[1] + alpha[1])
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
        oci2 <- setdiff(oci, c((1:N) * TT, (TT * (0:(N - 1)) + 1)))
        ## randomly select 1/3
        subcount <- floor(count / 3)
        rm.id <- sample(oci2, subcount, replace = FALSE)
        rm.id.upper <- rm.id - 1
        rm.id.lower <- rm.id + 1

        cv.id1 <- unique(c(rm.id, rm.id.upper, rm.id.lower))
        pos <- sapply(1:length(cv.id1), function(i) cv.id1[i] %in% oci)
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
# cv.sample <- function(I, D, count,
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
# }

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

    if (cv.count == 1 || count <= 2) { ## randomly missing
        cv.id <- sample(oci, count, replace = FALSE)
        rm.id.use <- cv.id
    } else {
        res <- TT %% cv.count
        int <- floor(TT / cv.count)
        rm.pos <- c()

        ## randomly select 1/cv.count
        subcount <- floor(count / cv.count * prop)

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
            for (i in 0:(cv.count - 1)) {
                rm.id.all <- c(rm.id.all, rm.id + i)
                if (i >= cv.donut && i <= (cv.count - cv.donut - 1)) {
                    rm.id.use <- c(rm.id.use, rm.id + i)
                }
            }
        }

        ## rm.id.all <- c(rm.id, rm.id + 1, rm.id + 2)

        rm.id.all <- intersect(rm.id.all, oci) ## remove missing values
        rm.id.all <- sort(rm.id.all)

        rm.id.use <- intersect(rm.id.use, oci)
        rm.id.use <- sort(rm.id.use)

        cv.id2 <- NULL
        if (length(rm.id.all) >= count) {
            cv.id <- rm.id.all[1:count]
            pos.cv <- which(rm.id.use >= min(cv.id) & rm.id.use <= max(cv.id))
            rm.id.use <- rm.id.use[pos.cv]
            # if (length(rm.id.use) >= count) {
            #    rm.id.use <- rm.id.use[1:count]
            # }
        } else {
            cv.id2 <- sample(setdiff(oci, rm.id.all), (count - length(rm.id.all)), replace = FALSE)
            cv.id <- sort(c(rm.id.all, cv.id2))
            ## rm.id.use <- sort(c(rm.id.use, cv.id2))
        }

        ## remove boundary observation
        # if (cv.treat == FALSE) {
        #    rm.pos <- c()
        #    for (i in 1:(cv.count - 1)) {
        #        rm.pos <- c(rm.pos, (TT * (0:(N-1)) + i))
        #    }
        #    oci2 <- setdiff(oci, rm.pos)
        # } else {
        #    rm.pos <- c()
        #    for (i in 1:(cv.count - 1)) {
        #        rm.pos <- c(rm.pos, (TT * (tr.pos - 1) + i))
        #    }
        #    oci2 <- setdiff(oci, rm.pos)
        # }

        ## randomly select 1/cv.count
        # subcount <- floor(count/cv.count)
        # if (subcount == 0) {
        #    subcount <- 1
        # }
        # rm.id <- sample(oci2, subcount, replace = FALSE)
        # rm.id.all <- c()
        # for (i in 1:(cv.count - 1)) {
        #    rm.id.all <- c(rm.id.all, rm.id - i)
        # }
        # rm.id.all <- c(rm.id.all, rm.id)

        # cv.id1 <- unique(rm.id.all)
        # pos <- unlist(sapply(1:length(cv.id1), function(i) cv.id1[i] %in% oci))
        # cv.id1 <- cv.id1[pos]

        # if (length(cv.id1) >= count) {
        #    cv.id1 <- cv.id1[1:count]
        #    cv.id2 <- NULL
        # } else {
        #    cv.id2 <- sample(setdiff(oci, cv.id1), (count - length(cv.id1)), replace = FALSE)
        # }

        # cv.id <- sort(c(cv.id1, cv.id2))
    }

    return(list(
        cv.id = cv.id, ## marked id
        est.id = rm.id.use
    )) ## id used for mspe
}

res.vcov <- function(res, ## TT*Nboots
                     cov.ar = 1) {
    T <- dim(res)[1]
    I <- is.na(res)
    count <- matrix(NA, T, T)

    res[is.na(res)] <- 0
    vcov <- res %*% t(res)


    for (i in 1:T) {
        for (j in 1:T) {
            if (i > j) {
                count[i, j] <- count[j, i]
            } else {
                if ((j - i) <= cov.ar) {
                    II <- I[i, ] + I[j, ]
                    count[i, j] <- min(1 / sum(II == 0), 1)
                } else {
                    count[i, j] <- 0
                }
            }
        }
    }
    vcov <- vcov * count
    return(vcov)
}


## ---------------------------------------------------------------------------
## Rolling (forward-only) cross-validation fold construction.
##
## For each unit (column of II), masks the LAST `cv.count` observed positions
## (in time order). Training uses everything else: the unit's earlier
## observations + all other units' full data. This closes the forward-leakage
## channel that the random-anchor cv.sample() leaves open at cv.donut = 0/1
## under serially correlated residuals.
##
## Arguments
##   II         TT x N integer matrix; II[t, i] = 1 if unit i is observed at t.
##   D          TT x N treatment indicator (0/1). When cv.treat = FALSE
##              (default for rolling on controls), D is ignored.
##   cv.count   integer; number of observations to mask per unit (default 3).
##   cv.treat   logical; if TRUE restrict masking to treated units only.
##              Default FALSE (mask control units).
##
## Returns a list with the same shape as cv.sample():
##   cv.id      integer vector of masked positions (1-indexed into vec(II))
##   rm.id.use  identical to cv.id (no donut shaving — every masked position
##              is scored, since rolling already prevents forward leakage).
##
## Deterministic: no random sampling. The same II + cv.count always produces
## the same cv.id. Callers therefore typically use k = 1 fold.
## ---------------------------------------------------------------------------
cv.sample.rolling <- function(II, D = NULL,
                              cv.count = 3,
                              cv.treat = FALSE,
                              min.T0 = 5L) {
    TT <- dim(II)[1]
    N  <- dim(II)[2]
    cv.id <- integer(0)

    ## Identify units to mask.
    if (isTRUE(cv.treat) && !is.null(D)) {
        ever_treat <- apply(D, 2, function(v) any(!is.na(v) & v >= 1))
        unit_idx   <- which(ever_treat)
    } else {
        unit_idx   <- seq_len(N)
    }

    ## Only mask units that retain >= min.T0 observations after masking.
    ## (Otherwise the downstream con2 check restores them and the rolling
    ## fold collapses to a tiny set of effective masks.)
    threshold <- min.T0 + cv.count
    for (j in unit_idx) {
        obs_t <- which(II[, j] == 1L)
        if (length(obs_t) >= threshold) {
            mask_t <- utils::tail(obs_t, cv.count)
            ## Convert (t, j) to a 1-indexed position into vec(II).
            cv.id <- c(cv.id, (j - 1L) * TT + mask_t)
        }
    }
    list(cv.id = cv.id, rm.id.use = cv.id)
}


v_replace <- function(needle, haystack) {
    sieved <- which(haystack == needle[1L])
    for (i in seq.int(1L, length(needle) - 1L)) {
        sieved <- sieved[haystack[sieved + i] == needle[i + 1L]]
    }
    out <- rep(NA, length(haystack))
    if (length(sieved) == 0) {
        return(out)
    }
    for (index in sieved) {
        if (!is.na(index)) {
            out[c(index:(index + length(needle) - 1))] <- needle
        }
    }
    return(out)
}


## Formula parser shared by fect.formula() and interFE.formula().
## fect reads a formula as a list of column names. Before 2.4.6 it used
## all.vars(), so a transformed term was silently fitted as its raw column:
## log(Y + 20) ~ D fitted Y ~ D, factor(g) entered as one linear slope,
## X1 * X2 as X1 + X2. Now every term must be a bare column name, and anything
## else stops with a message that says to create the column first.
## Intercept specifiers (+ 0, 0 +, + 1, 1 +, - 0, - 1, unary -1) are accepted
## and ignored, as before: fect always adds its own fixed effects.
## Returns list(Y = <outcome name>, rhs = <right-hand-side names, deduplicated,
## in order>) or stops.
.fect_formula_names <- function(formula, fun = "fect") {
    if (!inherits(formula, "formula") || length(formula) != 3L) {
        stop(fun, "() needs a two-sided formula such as Y ~ D + X1 + X2.",
             call. = FALSE)
    }
    is_icpt <- function(e) {
        is.numeric(e) && length(e) == 1L && e %in% c(0, 1)
    }
    flat <- function(e) {
        if (is.call(e) && identical(e[[1L]], as.name("+")) && length(e) == 3L) {
            return(c(flat(e[[2L]]), flat(e[[3L]])))
        }
        if (is.call(e) && identical(e[[1L]], as.name("-")) && length(e) == 3L &&
            is_icpt(e[[3L]])) {
            return(flat(e[[2L]]))                     # `... - 1`, `... - 0`
        }
        if (is.call(e) && identical(e[[1L]], as.name("-")) && length(e) == 2L &&
            is_icpt(e[[2L]])) {
            return(list())                            # unary `-1`, `-0`
        }
        if (is_icpt(e)) {
            return(list())                            # `+ 0`, `0 +`, `+ 1`, `1 +`
        }
        list(e)
    }
    is_col <- function(e) is.name(e) && !identical(as.character(e), ".")
    lhs <- formula[[2L]]
    rhs <- flat(formula[[3L]])
    bad <- c(
        if (!is_col(lhs)) deparse1(lhs),
        vapply(Filter(Negate(is_col), rhs), deparse1, "")
    )
    if (length(bad) > 0L) {
        stop(
            fun, "() formulas take bare column names only; ",
            if (length(bad) == 1L) "not a column name: " else "not column names: ",
            paste0("`", bad, "`", collapse = ", "), ". ",
            "Create the variable first (for example ",
            "`data$logY <- log(data$Y + 20)`) and use it in the formula. ",
            "For a categorical covariate create numeric dummy columns first ",
            "(for example with `model.matrix()`); for an interaction create ",
            "the product column first.",
            call. = FALSE
        )
    }
    if (length(rhs) == 0L) {
        stop(
            fun, "() needs ",
            if (identical(fun, "interFE")) {
                "at least one covariate on the right-hand side of the formula, for example Y ~ X1 + X2."
            } else {
                "a treatment variable on the right-hand side of the formula, for example Y ~ D + X1."
            },
            call. = FALSE
        )
    }
    list(Y = as.character(lhs),
         rhs = unique(vapply(rhs, as.character, "")))
}


## Covariates that cannot be estimated on the cells used to fit the model.
## Before 2.4.6 fect had no rank check: an exactly collinear covariate made
## X'X singular (garbage estimates, or "inv(): matrix is singular"), and a
## covariate absorbed by the fixed effects stopped with swapped labels.
## fect.default() calls this once on the final estimation panel and drops
## what it reports, as lm() does for aliased columns.
##
## X: TT x N x p covariate array. M: TT x N logical mask of the estimation
## cells. force: 0 none, 1 unit, 2 time, 3 two-way. extra.fe: TT x N x k
## array of additional fixed-effect codes (cfe), or NULL. Demeaning is
## unweighted: exact collinearity does not depend on positive weights.
## Returns list(drop = <indices into 1:p>, message = <one warning text or
## NULL>).
.fect_check_covariates <- function(X, M, force, Xname, extra.fe = NULL) {
    p <- dim(X)[3]
    none <- list(drop = integer(0), message = NULL)
    if (is.null(p) || p == 0) {
        return(none)
    }
    cells <- which(M)
    Xm <- matrix(NA_real_, length(cells), p)
    for (j in seq_len(p)) {
        Xm[, j] <- X[, , j][cells]
    }
    unit <- col(M)[cells]
    time <- row(M)[cells]
    fe <- list()
    if (force %in% c(1, 3)) fe$unit <- unit
    if (force %in% c(2, 3)) fe$time <- time
    if (!is.null(extra.fe) && dim(extra.fe)[3] > 0) {
        for (k in seq_len(dim(extra.fe)[3])) {
            fe[[paste0("extra", k)]] <- extra.fe[, , k][cells]
        }
    }
    ok <- stats::complete.cases(Xm)
    if (length(fe) > 0) {
        for (f in fe) ok <- ok & !is.na(f)
    }
    Xm <- Xm[ok, , drop = FALSE]
    fe <- lapply(fe, function(f) f[ok])
    unit <- unit[ok]
    time <- time[ok]
    if (nrow(Xm) == 0) {
        return(none)
    }

    reason <- rep(NA_character_, p)
    norm2 <- function(v) sqrt(sum(v^2))
    ## 1. no variation on the estimation cells (collinear with the intercept)
    c0 <- apply(Xm, 2, function(v) norm2(v - mean(v)))
    flat <- c0 <= 1e-12 * pmax(1, apply(Xm, 2, norm2))
    reason[flat] <- "has no variation on the cells used to estimate the covariate coefficients"

    ## 2. absorbed by the fixed effects
    rest <- which(!flat)
    if (length(rest) > 0) {
        if (length(fe) == 0) {
            Xd <- sweep(Xm, 2, colMeans(Xm))
        } else {
            Xd <- fixest::demean(Xm, f = fe, tol = 1e-10, iter = 10000,
                                 nthreads = 1, notes = FALSE)
            Xd <- matrix(Xd, nrow = nrow(Xm))
        }
        for (j in rest) {
            if (norm2(Xd[, j]) > 1e-8 * c0[j]) next
            if (force %in% c(1, 3) &&
                norm2(Xm[, j] - stats::ave(Xm[, j], unit)) <= 1e-8 * c0[j]) {
                reason[j] <- "does not vary over time within units, so it is absorbed by the unit fixed effects"
            } else if (force %in% c(2, 3) &&
                       norm2(Xm[, j] - stats::ave(Xm[, j], time)) <= 1e-8 * c0[j]) {
                reason[j] <- "does not vary across units within periods, so it is absorbed by the time fixed effects"
            } else if (force == 3 && is.null(extra.fe)) {
                reason[j] <- "is absorbed by the unit and time fixed effects together"
            } else {
                reason[j] <- "is absorbed by the fixed effects"
            }
        }
        ## 3. exact linear combinations of the other covariates; like lm(),
        ## keep the earlier column of a collinear set
        rest <- which(is.na(reason))
        if (length(rest) > 1) {
            q <- qr(Xd[, rest, drop = FALSE], tol = 1e-7)
            if (q$rank < length(rest)) {
                alias <- rest[q$pivot[-seq_len(q$rank)]]
                reason[alias] <- "is a linear combination of other covariates (after removing the fixed effects)"
            }
        }
    }

    drop <- which(!is.na(reason))
    if (length(drop) == 0) {
        return(none)
    }
    list(
        drop = drop,
        message = paste0(
            "Dropped ", length(drop), " covariate",
            if (length(drop) > 1) "s" else "",
            " that cannot be estimated on the cells used to fit the model: ",
            paste0("\"", Xname[drop], "\" ", reason[drop], collapse = "; "),
            ". ",
            if (length(drop) > 1) "Their coefficients are" else "Its coefficient is",
            " reported as NA."
        )
    )
}

## Put covariates dropped by .fect_check_covariates() back into a fit's
## outputs, so that every coefficient output has one row per requested
## covariate (NA for the dropped ones), and the stored covariate array and
## the internal estimator objects line up with the requested names again.
## keep: indices (into 1:p.all) of the covariates that were fitted, in order
## (empty when every covariate was dropped).
.fect_restore_covariates <- function(out, keep, p.all, X.all, se, binary) {
    ## p.all-row version of a length(keep)-row matrix `m`; all NA (with
    ## columns `cols`) when `m` does not have one row per kept covariate,
    ## e.g. the scalar NA the estimators return when no covariate is left
    expand <- function(m, cols) {
        m <- if (is.null(m)) NULL else as.matrix(m)
        if (length(keep) > 0 && !is.null(m) && nrow(m) == length(keep)) {
            full <- matrix(NA_real_, p.all, ncol(m),
                           dimnames = list(NULL, colnames(m)))
            full[keep, ] <- m
        } else {
            full <- matrix(NA_real_, p.all, length(cols),
                           dimnames = list(NULL, cols))
        }
        full
    }
    is.binary <- isTRUE(as.logical(binary))
    out$beta <- expand(out$beta, "Coef")
    if (is.binary) {
        out$marginal <- expand(out$marginal, "marginal")
    }
    if (isTRUE(as.logical(se))) {
        out$est.beta <- expand(
            out$est.beta, c("Coef", "S.E.", "CI.lower", "CI.upper", "p.value")
        )
        if (is.binary) {
            out$est.marginal <- expand(
                out$est.marginal,
                c("marginal", "S.E.", "CI.lower", "CI.upper", "p.value")
            )
        }
        if (!is.null(out$beta.boot) && length(keep) > 0) {
            bb <- as.matrix(out$beta.boot)
            if (nrow(bb) == length(keep)) {
                full <- matrix(NA_real_, p.all, ncol(bb))
                full[keep, ] <- bb
                out$beta.boot <- full
            }
        } else if (length(keep) == 0 && !is.null(out$att.avg.boot)) {
            out$beta.boot <- matrix(NA_real_, p.all, length(out$att.avg.boot))
        }
    }
    ## the estimator's covariate array and coefficient vectors, read by
    ## name-based lookups (plot(type = "hte"/"calendar"), fect_iden())
    if (is.array(out$X) && length(dim(out$X)) == 3 &&
        dim(out$X)[3] == length(keep)) {
        out$X <- X.all
    }
    for (slot in c("est", "est.cm")) {
        b <- out[[slot]]$beta
        if (length(keep) > 0 && !is.null(b) && length(b) == length(keep)) {
            full <- matrix(NA_real_, p.all, 1)
            full[keep, 1] <- as.numeric(b)
            out[[slot]]$beta <- full
        }
    }
    out
}
