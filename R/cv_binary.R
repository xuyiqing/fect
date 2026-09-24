## Binary IFE: explicit training masks and probability scoring.

## Bai-Ng style information criterion for the Probit IFE model, as reported
## by fect 1.1.x ("IC" in its binary CV table):
##   PC(r) = r * (N + T) / (N T) * log(N T / (N + T)) - 2 * mean log-likelihood,
## where the mean log-likelihood is averaged over the cells used in the fit.
## It is a function of the full-sample fit only, so it does not depend on
## how the data are split for cross-validation. The C++ `IC` slot is the
## BIC-type penalty and is kept separately.
.fect_binary_pc <- function(loglikelihood, r, N, TT) {
    r * (N + TT) / (N * TT) * log(N * TT / (N + TT)) - 2 * loglikelihood
}
.fect_binary_estimate <- function(Y, X, II, r, force, QR, tol,
                                  max.iteration = 1000) {
    TT <- nrow(Y); N <- ncol(Y)
    p <- if (is.null(X)) 0L else dim(X)[3L]
    data.ini <- cbind(c(Y), rep(seq_len(N), each=TT), rep(seq_len(TT), N))
    if (p > 0L) data.ini <- cbind(data.ini, matrix(X, TT*N, p))
    oci <- which(II == 1L)
    if (length(unique(Y[oci])) != 2L)
        stop("Binary training outcomes must contain both 0 and 1.")
    # Excluded outcomes cannot enter initialization, including PCA residuals.
    data.ini[-oci, 1L] <- 0
    init <- BiInitialFit(data.ini, QR=QR, r=r, force=force, oci=oci)
    YY <- Y; YY[II == 0L] <- 0
    XX <- if (p > 0L) X else array(0, c(TT,N,0))
    if (isTRUE(as.logical(QR))) {
        est <- inter_fe_d_qr_ub(YY, init$Y0, init$FE0, init$factor0,
            init$xi0, XX, II, r, force, tol=tol, mniter=max.iteration)
    } else {
        est <- inter_fe_d_ub(YY, init$Y0, init$FE0, XX, II, r, force,
                            tol=tol, mniter=max.iteration)
    }
    if (any(!is.finite(est$fit))) stop("Binary fit returned nonfinite indices.")
    est
}

.fect_binary_folds <- function(II, D, k, cv.method, cv.nobs, cv.buffer,
                               cv.prop, min.T0, cv.donut, seed=NULL) {
    if (cv.method == "rolling")
        return(.build_cv_mask_rolling(II, D, k, cv.nobs, cv.buffer, cv.prop, min.T0, seed=seed))
    count <- max(1L, floor(sum(II)*cv.prop))
    folds <- vector("list", k)
    for (i in seq_len(k)) {
        valid <- FALSE
        for (attempt in seq_len(100L)) {
            f <- cv.sample(II, D, count, cv.nobs,
                cv.treat=cv.method == "treated_units", cv.donut=cv.donut)
            train <- II; train[f$cv.id] <- 0L
            valid <- length(f$est.id) > 0L &&
                all(colSums(train) >= min.T0) && all(rowSums(train) > 0L)
            if (valid) break
        }
        if (!valid) stop("Cannot construct binary CV folds with the requested training history and donor support.")
        folds[[i]] <- f
    }
    folds
}

fect_binary_cv <- function(Y, X, D, I, II, T.on, T.off=NULL,
    k=20, cv.prop=0.1, cv.method="rolling", cv.nobs=3, cv.buffer=1,
    min.T0=5, r=0, r.end, QR=FALSE, force, hasRevs=0, tol,
    group.level=NULL, group=NULL, cv.donut=1, cv.rule="1se", max.iteration=1000, seed=NULL) {
    cv.method <- .fect_normalize_cv_method(cv.method,
        allowed=c("rolling","block","all_units","treated_units"))
    cv.rule <- .fect_validate_cv_rule(cv.rule)
    if (length(k) != 1L || !is.finite(k) || k < 1 || k != as.integer(k))
        stop("k must be a positive integer.")
    if (!is.finite(cv.prop) || cv.prop <= 0 || cv.prop > 1)
        stop("cv.prop must satisfy 0 < cv.prop <= 1.")
    if (r.end < r || r < 0) stop("Invalid binary candidate rank range.")
    if (!is.null(seed)) set.seed(seed)
    ranks <- seq.int(r, r.end)
    folds <- .fect_binary_folds(II,D,k,cv.method,cv.nobs,cv.buffer,cv.prop,min.T0,cv.donut,seed)
    losses <- classification <- matrix(NA_real_, length(ranks), k)
    failures <- matrix(NA_character_, length(ranks), k)
    counts <- vapply(folds, function(f) length(f$est.id), integer(1))
    full <- vector("list", length(ranks))
    full.failures <- rep(NA_character_, length(ranks))
    p <- if (is.null(X)) 0L else dim(X)[3L]
    feasible <- function(mask, rank) {
        all(colSums(mask) >= max(min.T0,rank+1L)) &&
        all(rowSums(mask) >= rank+1L) &&
        sum(mask) > rank*(nrow(mask)+ncol(mask)-rank)+p
    }
    message("Cross-validating binary probability predictions ...")
    for (a in seq_along(ranks)) {
        rank <- ranks[a]
        for (b in seq_len(k)) {
            f <- folds[[b]]; mask <- II; mask[f$cv.id] <- 0L
            fit <- tryCatch({
                if (!feasible(mask,rank)) stop("insufficient training history or donor support for rank")
                z <- .fect_binary_estimate(Y,X,mask,rank,force,QR,tol,max.iteration)
                if (z$niter >= max.iteration) stop("iteration limit reached")
                z
            }, error=function(e)e)
            if (inherits(fit,"error")) {
                failures[a,b] <- conditionMessage(fit); next
            }
            predicted <- pnorm(fit$fit[f$est.id])
            losses[a,b] <- mean((Y[f$est.id]-predicted)^2)
            classification[a,b] <- mean(Y[f$est.id] != (predicted >= 0.5))
        }
        full[[a]] <- tryCatch({
            if (!feasible(II,rank)) stop("insufficient full-sample support for rank")
            z <- .fect_binary_estimate(Y,X,II,rank,force,QR,tol,max.iteration)
            if (z$niter >= max.iteration) stop("iteration limit reached")
            z
        }, error=function(e)e)
        if (inherits(full[[a]],"error")) full.failures[a] <- conditionMessage(full[[a]])
    }
    # Compare candidates on every common fold; failures cannot improve scores.
    eligible <- rowSums(is.finite(losses)) == k & is.na(full.failures)
    means <- rowMeans(losses)
    ses <- if (k > 1L) apply(losses,1L,stats::sd)/sqrt(k) else rep(0,length(ranks))
    means[!eligible] <- Inf
    chosen <- .fect_apply_cv_rule(means,ses,rule=cv.rule)
    if (is.na(chosen)) stop("All binary CV candidates failed: ",
        paste(unique(na.omit(c(failures,full.failures))),collapse="; "))
    if (any(!eligible)) warning("Some binary CV ranks failed or lacked support; see cv.failures.",call.=FALSE)
    # Share the current result builder with fixed-rank and inference paths.
    out <- fect_fe(Y=Y,X=X,D=D,W=NULL,I=I,II=II,T.on=T.on,T.off=T.off,
        r.cv=ranks[chosen],binary=TRUE,QR=QR,force=force,hasRevs=hasRevs,tol=tol,
        max.iteration=max.iteration,group.level=group.level,group=group)
    get_full <- function(name) vapply(full,function(z) {
        if (inherits(z,"error")) NA_real_ else as.numeric(z[[name]])
    },numeric(1))
    loglik.full <- get_full("loglikelihood")
    PC <- .fect_binary_pc(loglik.full, ranks, ncol(Y), nrow(Y))
    out$CV.out <- cbind(r=ranks,IC=get_full("IC"),PC=PC,
        `Log-likelihood`=loglik.full,MSPE=means,
        MSPE.SE=ses,Classification.error=rowMeans(classification))
    for (a in seq_along(ranks)) {
        message(sprintf(" r = %d; PC = %.5f; Log-likelihood = %.5f; MSPE = %.5f%s",
            ranks[a], PC[a], loglik.full[a], means[a],
            if (a == chosen) " *" else ""))
    }
    out$cv.rule <- cv.rule
    out$cv.method <- if (cv.method == "all_units") "block" else cv.method
    out$cv.loss <- "probability_mspe"
    out$cv.folds <- folds
    out$cv.loss.per.fold <- losses
    out$cv.classification.per.fold <- classification
    out$cv.counts <- counts
    out$cv.pooled.mspe <- as.numeric(losses %*% counts/sum(counts))
    out$cv.failures <- list(folds=failures,full=full.failures)
    out$cv.settings <- list(k=k,cv.prop=cv.prop,cv.nobs=cv.nobs,
        cv.buffer=cv.buffer,min.T0=min.T0,cv.donut=cv.donut)
    # Units the folds could draw from, as the fold builder counted them: only
    # untreated cells before the first treated period count, so cells after
    # a reversal do not make a unit eligible.
    if (cv.method == "rolling") out$cv.eligible.units <- attr(folds, "eligible.units")
    message("Selected r = ",out$r.cv," (",cv.rule,", probability MSPE).")
    out
}
