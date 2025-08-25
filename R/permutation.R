###################################################################
## permutation
###################################################################
fect_permu <- function(Y,
                       X,
                       D,
                       I,
                       permu.dimension = "time",
                       r.cv = NULL,
                       lambda.cv = NULL,
                       m = 2,
                       method = "ife",
                       degree = 2,
                       knots = NULL,
                       force,
                       tol,
                       norm.para,
                       nboots,
                       parallel = TRUE,
                       cores = NULL) {
    TT <- dim(Y)[1]
    n.shuffle <- TT %/% m
    t.pos <- rep(1:n.shuffle, each = m)
    if (TT > n.shuffle * m) {
        t.pos <- c(t.pos, rep(n.shuffle + 1, TT - n.shuffle * m))
    }
    tt <- 1:TT
    ## list
    l.tt <- split(tt, t.pos)

    if (permu.dimension == "time") {
        sub.permu <- function() {
            sub.pos <- c()
            tt.length <- length(l.tt)
            one.rank <- sample(1:tt.length, tt.length, replace = FALSE)
            for (i in 1:tt.length) {
                sub.pos <- c(sub.pos, unlist(l.tt[[one.rank[i]]]))
            }
            Y.permu <- as.matrix(Y[sub.pos, ])
            I.permu <- as.matrix(I[sub.pos, ])
            # D.permu <- as.matrix(D[sub.pos, ])
            if (!is.null(X)) {
                X.permu <- X[sub.pos, , , drop = FALSE]
            }
            result <- try(one.permu(
                Y.permu,
                X.permu,
                D,
                I.permu,
                r.cv,
                lambda.cv,
                method,
                degree,
                knots,
                force,
                tol,
                norm.para
            ), silent = TRUE)


            if ("try-error" %in% class(result)) {
                return(NA)
            } else {
                return(result)
            }
        }
    }

    if (permu.dimension == "unit") {
        NN <- dim(Y)[2]
        index.NN <- c(1:NN)

        sub.permu <- function() {
            sample.unit <- sample(index.NN,
                size = length(index.NN),
                replace = FALSE
            )

            Y.permu <- as.matrix(Y[, sample.unit])
            I.permu <- as.matrix(I[, sample.unit])
            if (!is.null(X)) {
                X.permu <- X[, sample.unit, , drop = FALSE]
            }
            result <- try(one.permu(
                Y.permu,
                X.permu,
                D,
                I.permu,
                r.cv,
                lambda.cv,
                method,
                degree,
                knots,
                force,
                tol,
                norm.para
            ), silent = TRUE)


            if ("try-error" %in% class(result)) {
                return(NA)
            } else {
                return(result)
            }
        }
    }


    att.avg.permu <- rep(NA, nboots)

    if (parallel == TRUE) {
        permu.out <- foreach(
            j = 1:nboots,
            .inorder = FALSE,
            .export = c("one.permu"),
            .packages = c("fect")
        ) %dopar% {
            return(sub.permu())
        }

        for (j in 1:nboots) {
            att.avg.permu[j] <- permu.out[[j]]
        }
    } else {
        for (j in 1:nboots) {
            permu <- sub.permu()
            att.avg.permu[j] <- permu
            ## report progress
            if (j %% 100 == 0) {
                message(".")
            }
        }
    }

    if (sum(is.na(c(att.avg.permu))) > 0) {
        permu.rm <- which(is.na(c(att.avg.permu)))
        att.avg.permu <- att.avg.permu[-permu.rm]
    }

    message(length(att.avg.permu), " permutes\n", sep = "")

    return(att.avg.permu)
}


## sub-function

one.permu <- function(Y, # Outcome variable, (T*N) matrix
                      X, # Explanatory variables:  (T*N*p) array
                      D, #  Indicator for treated unit (tr==1)
                      I,
                      r.cv = 0, # initial number of factors considered if CV==1
                      lambda.cv = 1,
                      method = "fe",
                      degree = 2,
                      knots = NULL,
                      force,
                      tol, # tolerance level
                      norm.para = NULL) {
    ## -------------------------------##
    ## Parsing data
    ## -------------------------------##
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

    ## observed indicators:
    II <- I
    II[which(D == 1)] <- 0

    ## remove period when no obserbations are under control
    T0 <- apply(II, 1, sum)
    XX <- NULL
    if (min(T0) == 0) {
        t0.pos <- which(T0 == 0)
        TT <- TT - length(t0.pos)
        if (TT <= 2) {
            stop("No enough observations under control.\n")
        }
        Y <- as.matrix(Y[-t0.pos, ])
        D <- as.matrix(D[-t0.pos, ])
        I <- as.matrix(I[-t0.pos, ])
        II <- as.matrix(II[-t0.pos, ])
        if (is.null(X) == FALSE) {
            XX <- array(NA, dim = c(TT, N, p))
            for (i in 1:p) {
                subX <- X[, , i]
                XX[, , i] <- as.matrix(subX[-t0.pos, ])
            }
            X <- XX
        }
    }

    ## replicate data
    YY <- Y
    YY[which(II == 0)] <- 0 ## reset to 0
    D[which(I == 0)] <- 0
    oci <- which(c(II) == 1)

    if (sum(D) == 0) {
        stop("No valid observations under treatment.\n")
    }

    if (method %in% c("polynomial", "bspline")) {
        stop("Doesn't support in this version.")
    } else {
        ## initial fit using fastplm
        data.ini <- matrix(NA, (TT * N), (2 + 1 + p))
        data.ini[, 2] <- rep(1:N, each = TT) ## unit fe
        data.ini[, 3] <- rep(1:TT, N) ## time fe
        data.ini[, 1] <- c(Y) ## outcome
        if (p > 0) { ## covar
            for (i in 1:p) {
                data.ini[, (3 + i)] <- c(X[, , i])
            }
        }
        ## observed Y0 indicator:
        initialOut <- Y0 <- beta0 <- FE0 <- xi0 <- factor0 <- NULL

        initialOut <- initialFit(data = data.ini, force = force, oci = oci)
        Y0 <- initialOut$Y0
        beta0 <- initialOut$beta0
        if (p > 0 && sum(is.na(beta0)) > 0) {
            beta0[which(is.na(beta0))] <- 0
        }

        ## -------------------------------##
        ## ----------- Main Algorithm ----------- ##
        ## -------------------------------##

        est <- NULL
        if (method == "fe") {
            est <- inter_fe_ub(YY, Y0, X, II, beta0, 0, force = force, tol)
        } else if (method == "ife") {
            est <- inter_fe_ub(YY, Y0, X, II, beta0, r.cv, force = force, tol)
        } else if (method == "mc") {
            est <- inter_fe_mc(YY, Y0, X, II, beta0, 1, lambda.cv, force, tol)
        }

        if (!is.null(norm.para)) {
            est$fit <- est$fit * norm.para[1]
        }
        Y.ct <- est$fit
    }

    if (!is.null(norm.para)) {
        Y <- Y * norm.para[1]
    }

    eff <- Y - Y.ct

    att.avg <- sum(eff * D) / (sum(D))

    return(abs(att.avg))
}
