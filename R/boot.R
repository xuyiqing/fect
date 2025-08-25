###############################################
## Inference
###############################################
basic_ci_alpha <- function(theta, boots, alpha) {
  n_p <- length(theta)
  probs <- c(alpha/2, 1 - alpha/2)
  ci_mat <- t(sapply(seq_len(n_p), function(i) {
    vals <- boots[i, ]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 1) {
      qs <- tryCatch(
        quantile(vals, probs = probs, na.rm = TRUE),
        error = function(e) c(NA_real_, NA_real_)
      )
      if (all(!is.na(qs))) {
        # reflect about theta[i]
        lower <- 2 * theta[i] - qs[2]
        upper <- 2 * theta[i] - qs[1]
        return(c(lower, upper))
      }
    }
    c(NA_real_, NA_real_)
  }))

  colnames(ci_mat) <- c("lower", "upper")
  ci_mat
}

fect_boot <- function(Y,
                      X,
                      D,
                      W,
                      cl = NULL,
                      I,
                      II,
                      T.on,
                      T.off = NULL,
                      T.on.carry = NULL,
                      T.on.balance = NULL,
                      balance.period = NULL,
                      method = "ife",
                      degree = 2,
                      sfe = NULL,
                      cfe = NULL,
                      ind.matrix = NULL,
                      knots = NULL,
                      criterion = "mspe",
                      CV = 0,
                      k = 10,
                      cv.prop = 0.1,
                      cv.treat = TRUE,
                      cv.nobs = 3,
                      r = 0,
                      r.end = 3,
                      lambda = NULL,
                      nlambda = 10,
                      alpha = 0.05,
                      binary = 0,
                      QR = 0,
                      force = 0,
                      hasRevs = 1,
                      tol = 1e-3,
                      max.iteration = 1000,
                      norm.para = NULL,
                      placebo.period = NULL,
                      placeboTest = 0,
                      carryoverTest = 0,
                      carryover.period = NULL,
                      vartype = "bootstrap",
                      quantile.CI = FALSE,
                      nboots = 200,
                      parallel = TRUE,
                      cores = NULL,
                      group.level = NULL,
                      group = NULL,
                      dis = 0,
                      keep.sims = FALSE) {
  na.pos <- NULL
  TT <- dim(Y)[1]
  N <- dim(Y)[2]
  if (is.null(X) == FALSE) {
    p <- dim(X)[3]
  } else {
    p <- 0
  }
  if (is.null(W)) {
    W.use <- as.matrix(0)
    use_weight <- 0
  } else {
    W.use <- W
    W.use[which(II == 0)] <- 0
    use_weight <- 1
  }



  if (hasRevs == 1) {
    ## D.fake : check reversals
    D.fake <- apply(D, 2, function(vec) {
      cumsum(vec)
    })
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
    if (method == "gsynth") {
      out <- fect_gsynth(
        Y = Y, X = X, D = D, W = W, I = I, II = II,
        T.on = T.on, T.off = T.off, CV = 0,
        T.on.balance = T.on.balance,
        balance.period = balance.period,
        r = r, binary = binary, QR = QR,
        force = force, hasRevs = hasRevs,
        tol = tol, max.iteration = max.iteration, boot = 0,
        norm.para = norm.para,
        placebo.period = placebo.period,
        placeboTest = placeboTest,
        carryover.period = carryover.period,
        carryoverTest = carryoverTest,
        group.level = group.level, group = group
      )
    } else if (method == "ife") {
      out <- fect_fe(
        Y = Y, X = X, D = D, W = W, I = I, II = II,
        T.on = T.on, T.off = T.off, T.on.carry = T.on.carry,
        T.on.balance = T.on.balance,
        balance.period = balance.period,
        r.cv = r, binary = binary, QR = QR,
        force = force, hasRevs = hasRevs,
        tol = tol, max.iteration = max.iteration, boot = 0,
        norm.para = norm.para,
        placebo.period = placebo.period,
        placeboTest = placeboTest,
        carryover.period = carryover.period,
        carryoverTest = carryoverTest,
        group.level = group.level, group = group
      )
    } else if (method == "mc") {
      out <- try(fect_mc(
        Y = Y, X = X, D = D, W = W, I = I, II = II,
        T.on = T.on, T.off = T.off, T.on.carry = T.on.carry,
        T.on.balance = T.on.balance,
        balance.period = balance.period,
        lambda.cv = lambda, force = force, hasRevs = hasRevs,
        tol = tol, max.iteration = max.iteration, boot = 0,
        norm.para = norm.para,
        placebo.period = placebo.period,
        placeboTest = placeboTest,
        carryover.period = carryover.period,
        carryoverTest = carryoverTest,
        group.level = group.level, group = group
      ), silent = TRUE)
      if ("try-error" %in% class(out)) {
        stop("\nCannot estimate using full data with MC algorithm.\n")
      }
    } else if (method %in% c("polynomial", "bspline", "cfe")) {
      out <- try(fect_polynomial(
        Y = Y, D = D, X = X, W = W, I = I,
        II = II, T.on = T.on, T.on.carry = T.on.carry,
        T.on.balance = T.on.balance,
        balance.period = balance.period,
        T.off = T.off,
        method = method, degree = degree,
        knots = knots, force = force,
        sfe = sfe, cfe = cfe,
        ind.matrix = ind.matrix,
        hasRevs = hasRevs,
        tol = tol, max.iteration = max.iteration, boot = 0,
        placeboTest = placeboTest,
        placebo.period = placebo.period,
        carryover.period = carryover.period,
        carryoverTest = carryoverTest,
        norm.para = norm.para,
        group.level = group.level,
        group = group
      ), silent = TRUE)
      # I.report <- out$I
      # II.report <- out$II
      if ("try-error" %in% class(out)) {
        stop("\nCannot estimate.\n")
      }
    }
  } else {
    ## cross-valiadtion
    if (binary == 0) {
      out <- fect_cv(
        Y = Y, X = X, D = D, W = W, I = I, II = II,
        T.on = T.on, T.off = T.off, T.on.carry = T.on.carry,
        T.on.balance = T.on.balance,
        balance.period = balance.period,
        method = method, criterion = criterion,
        k = k, r = r, r.end = r.end,
        nlambda = nlambda, lambda = lambda,
        force = force, hasRevs = hasRevs,
        tol = tol, max.iteration = max.iteration, norm.para = norm.para,
        group.level = group.level, group = group,
        cv.prop = cv.prop, cv.treat = cv.treat,
        cv.nobs = cv.nobs
      )

      # method <- out$method
    } else {
      out <- fect_binary_cv(
        Y = Y, X = X, D = D,
        I = I, II = II,
        T.on = T.on, T.off = T.off,
        k = k, r = r, r.end = r.end,
        QR = QR, force = force,
        hasRevs = hasRevs, tol = tol,
        group.level = group.level, group = group
      )
      # method <- "ife"
    }
  }



  ## output
  validX <- out$validX
  eff <- out$eff
  att.avg <- out$att.avg
  att.avg.unit <- out$att.avg.unit
  calendar.eff <- out$eff.calendar
  calendar.eff.fit <- out$eff.calendar.fit
  calendar.N <- out$N.calendar

  group.att <- out$group.att

  att <- out$att
  time.on <- out$time
  target.enp <- out$calendar.enp

  time.off <- NULL
  if (hasRevs == 1) {
    att.off <- out$att.off
    time.off <- out$time.off
  }
  carry.att <- carry.time <- NULL
  if (!is.null(T.on.carry)) {
    carry.att <- out$carry.att
    carry.time <- out$carry.time
  }

  if (!is.null(balance.period)) {
    balance.att <- out$balance.att
    balance.time <- out$balance.time
    balance.count <- out$balance.count
    balance.avg.att <- out$balance.avg.att
    if (!is.null(placebo.period) & placeboTest == TRUE) {
      balance.att.placebo <- out$balance.att.placebo
    }
  }

  if (!is.null(W)) {
    att.avg.W <- out$att.avg.W
    att.on.sum.W <- out$att.on.sum.W
    att.on.W <- out$att.on.W
    count.on.W <- out$count.on.W
    time.on.W <- out$time.on.W
    W.on.sum <- out$W.on.sum
    if (!is.null(placebo.period) & placeboTest == TRUE) {
      att.placebo.W <- out$att.placebo.W
    }
    if (hasRevs == 1) {
      att.off.sum.W <- out$att.off.sum.W
      att.off.W <- out$att.off.W
      count.off.W <- out$count.off.W
      time.off.W <- out$time.off.W
      W.off.sum <- out$W.off.sum
    }
    if (!is.null(carryover.period) & carryoverTest == TRUE) {
      att.carryover.W <- out$att.carryover.W
    }
  }

  eff.out <- out$eff
  fit.out <- out$Y.ct
  N_unit <- dim(out$res)[2]

  if (!is.null(group)) {
    group.output.origin <- out$group.output
    group.output.name <- names(out$group.output)
    group.time.on <- list()
    group.time.off <- list()
    for (sub.name in group.output.name) {
      group.time.on[[sub.name]] <- out$group.output[[sub.name]]$time.on
      if (hasRevs == 1) {
        group.time.off[[sub.name]] <- out$group.output[[sub.name]]$time.off
      }
    }
  }

  if (p > 0) {
    beta <- out$beta
  } else {
    beta <- matrix(0, 1, 0)
  }

  if (is.null(cl)) {
    cl.unique <- NULL
  } else {
    cl.unique <- unique(c(cl))
  }

  if (vartype == "jackknife") {
    nboots <- N
  }

  ## bootstrapped estimates
  if (keep.sims) {
    if (vartype=="jackknife") {
      D.boot <- array(NA, dim = c(TT, N-1, nboots))
      I.boot <- array(NA, dim = c(TT, N-1, nboots))
      eff.boot <- array(0, dim = c(TT, N-1, nboots))  ## to store results
    } else {
      D.boot <- array(NA, dim = c(TT, N, nboots))
      I.boot <- array(NA, dim = c(TT, N, nboots))
      eff.boot <- array(0, dim = c(TT, N, nboots))  ## to store results
    }
    colnames.boot <- c()
  }

  att.avg.boot <- matrix(0, 1, nboots)
  att.avg.unit.boot <- matrix(0, 1, nboots)
  att.boot <- matrix(0, length(time.on), nboots)
  att.count.boot <- matrix(0, length(time.on), nboots)
  beta.boot <- marginal.boot <- att.off.boot <- att.off.count.boot <- NULL
  calendar.eff.boot <- matrix(0, TT, nboots)
  calendar.eff.fit.boot <- matrix(0, TT, nboots)

  if (hasRevs == 1) {
    eff.off.boot <- c()
    att.off.boot <- matrix(0, length(time.off), nboots)
    att.off.count.boot <- matrix(0, length(time.off), nboots)
  }
  if (!is.null(T.on.carry)) {
    carry.att.boot <- matrix(0, length(carry.att), nboots)
  }
  if (!is.null(balance.period)) {
    balance.att.boot <- matrix(0, length(balance.att), nboots) # dynamic att
    balance.count.boot <- matrix(0, length(balance.att), nboots)
    balance.avg.att.boot <- matrix(0, 1, nboots)
    if (!is.null(placebo.period) & placeboTest == TRUE) {
      balance.att.placebo.boot <- matrix(0, 1, nboots)
    }
  }

  if (!is.null(W)) {
    att.avg.W.boot <- matrix(0, 1, nboots)
    att.on.W.boot <- matrix(0, length(time.on.W), nboots)
    att.on.count.W.boot <- matrix(0, length(time.on.W), nboots)
    if (!is.null(placebo.period) & placeboTest == TRUE) {
      att.placebo.W.boot <- matrix(0, 1, nboots)
    }
    if (hasRevs == 1) {
      att.off.W.boot <- matrix(0, length(time.off.W), nboots)
      att.off.count.W.boot <- matrix(0, length(time.off.W), nboots)
      if (!is.null(carryover.period) & carryoverTest == TRUE) {
        att.carryover.W.boot <- matrix(0, 1, nboots)
      }
    }
  }


  if (p > 0) {
    beta.boot <- matrix(0, p, nboots)
    if (binary == TRUE) {
      marginal.boot <- matrix(0, p, nboots)
    }
  }
  if (!is.null(placebo.period) & placeboTest == TRUE) {
    att.placebo.boot <- matrix(0, 1, nboots)
  }
  if (!is.null(carryover.period) & carryoverTest == TRUE) {
    att.carryover.boot <- matrix(0, 1, nboots)
  }

  group.att.boot <- NULL
  if (!is.null(group.att)) {
    group.att.boot <- matrix(0, length(group.att), nboots)
  }

  group.atts.boot <- NULL
  group.atts.off.boot <- NULL
  group.atts.placebo.boot <- NULL
  group.atts.carryover.boot <- NULL
  if (!is.null(group)) {
    group.atts.boot <- list()
    group.atts.off.boot <- list()
    group.att.placebo.boot <- list()
    group.att.carryover.boot <- list()
    for (sub.name in group.output.name) {
      subgroup.time.on <- group.time.on[[sub.name]]
      group.atts.boot[[sub.name]] <- matrix(0, length(subgroup.time.on), nboots)
      if (hasRevs == 1) {
        subgroup.time.off <- group.time.off[[sub.name]]
        group.atts.off.boot[[sub.name]] <- matrix(0, length(subgroup.time.off), nboots)
      }
      if (placeboTest) {
        group.att.placebo.boot[[sub.name]] <- matrix(0, 1, nboots)
      }
      if (carryoverTest) {
        group.att.carryover.boot[[sub.name]] <- matrix(0, 1, nboots)
      }
    }
  }

  if (dis) {
    if (vartype == "jackknife") {
      message("Jackknife estimates ... ")
    } else {
      message("Bootstrapping for uncertainties ... ")
    }
  }

  if (binary == TRUE & vartype == "parametric") {
    one.nonpara <- function(num = NULL) {
      Y.boot <- Y
      Y.fit <- out$Y.ct
      Y.fit[which(is.na(Y.fit))] <- 0
      Y.boot.new <- matrix(sapply(1:length(c(Y.fit)), function(i) {
        rbinom(1, 1, c(Y.fit)[i])
      }), TT, N)
      Y.boot[which(II == 1)] <- Y.boot.new[which(II == 1)]

      placebo.period.boot <- NULL
      carryover.period.boot <- NULL
      if (placeboTest == TRUE) {
        placebo.period.boot <- placebo.period
      }
      if (carryoverTest == TRUE) {
        carryover.period.boot <- carryover.period
      }

      boot <- try(fect_fe(
        Y = Y.boot, X = X, D = D, W = W,
        I = I, II = II,
        T.on = T.on, T.off = T.off, T.on.carry = T.on.carry,
        T.on.balance = T.on.balance,
        balance.period = balance.period,
        r.cv = out$r.cv, binary = binary,
        QR = QR, force = force,
        hasRevs = hasRevs, tol = tol, max.iteration = max.iteration, boot = 1,
        norm.para = norm.para,
        calendar.enp.seq = target.enp,
        time.on.seq = time.on,
        time.off.seq = time.off,
        time.on.carry.seq = carry.time,
        time.on.balance.seq = balance.time,
        time.on.seq.W = time.on.W,
        time.off.seq.W = time.off.W,
        placebo.period = placebo.period.boot,
        placeboTest = placeboTest,
        carryoverTest = carryoverTest,
        carryover.period = carryover.period.boot,
        group.level = group.level,
        group = group
      ), silent = TRUE)

      if ("try-error" %in% class(boot)) {
        boot0 <- list(
          att.avg = NA, att = NA, count = NA,
          beta = NA, att.off = NA, count.off = NA, eff.calendar = NA,
          eff.calendar.fit = NA,
          att.placebo = NA, att.avg.unit = NA, att.carryover = NA,
          group.att = NA, marginal = NA, carry.att = NA, balance.att = NA,
          balance.att.placebo = NA, balance.count = NA,
          att.avg.W = NA, att.on.W = NA, count.on.W = NA, time.on.W = NA, att.placebo.W = NA,
          att.off.W = NA, count.off.W = NA, time.off.W = NA, att.carryover.W = NA,
          balance.avg.att = NA, balance.time = NA, group.output = list()
        )
        return(boot0)
      } else {
        return(boot)
      }
    }
  } else if (binary == FALSE & method %in% c("gsynth") & vartype == "parametric") {
    message("Parametric Bootstrap \n")
    sum.D <- colSums(out$D)
    id.tr <- which(sum.D > 0)
    I.tr <- as.matrix(out$I[, id.tr])
    id.co <- which(sum.D == 0)
    Nco <- length(id.co)
    Ntr <- length(id.tr)

    fit.out[which(out$I == 0)] <- 0

    error.co <- out$res.full[, id.co]
    I.co <- out$I[, id.co]
    T0.ub <- apply(as.matrix(out$D[, id.tr] == 0), 2, sum)
    T0.ub.min <- min(T0.ub)
    co.pre <- apply(as.matrix(I.co[1:T0.ub.min, ]), 2, sum)
    co.post <- apply(as.matrix(I.co[(max(T0.ub) + 1):TT, ]), 2, sum)
    if (force %in% c(1, 3)) {
      valid.co <- id.co[(co.pre >= (out$r.cv + 1)) & (co.post >= 1)]
    } else {
      valid.co <- id.co[(co.pre >= out$r.cv) & (co.post >= 1)]
    }

    draw.error <- function() {
      repeat {
        fake.tr <- sample(id.co, 1, replace = FALSE)
        if (fake.tr %in% valid.co) {
          break
        }
      }

      id.co.rest <- id.co[which(!id.co %in% fake.tr)]
      repeat {
        id.co.pseudo <- sample(id.co.rest, Nco, replace = TRUE)
        if (sum(apply(as.matrix(out$I[, id.co.pseudo]), 1, sum) >= 1) == TT) {
          break
        }
      }

      id.pseudo <- c(rep(fake.tr, Ntr), id.co.pseudo) ## Ntr + ...
      I.id.pseudo <- out$I[, id.pseudo]
      II.id.pseudo <- out$II[, id.pseudo]
      ## obtain the prediction eror
      D.pseudo <- out$D[, c(id.tr, id.co.pseudo)] ## fake.tr + control left
      Y.pseudo <- out$Y[, id.pseudo]
      T.on.pseudo <- T.on[, id.pseudo]

      X.pseudo <- NULL
      if (p > 0) {
        X.pseudo <- X[, id.pseudo, , drop = FALSE]
      }

      ## output
      # synth.out <- try(fect_gsynth(Y = Y.pseudo, X = X.pseudo, D = D.pseudo,
      #                             I = I.id.pseudo, II = II.id.pseudo,
      #                             force = force, r = out$r.cv, CV = 0,
      #                             tol = tol, norm.para = norm.para, boot = 1), silent = TRUE)

      synth.out <- try(fect_gsynth(
        Y = Y.pseudo, X = X.pseudo, D = D.pseudo, W = NULL,
        I = I.id.pseudo, II = II.id.pseudo,
        T.on = T.on.pseudo, hasRevs = hasRevs,
        force = force, r = out$r.cv, CV = 0,
        tol = tol, max.iteration = max.iteration, norm.para = norm.para, boot = 1
      ), silent = TRUE)

      if ("try-error" %in% class(synth.out)) {
        return(matrix(NA, TT, Ntr))
      } else {
        if ("eff" %in% names(synth.out)) {
          if (is.null(norm.para)) {
            output <- synth.out$eff.tr
          } else {
            output <- synth.out$eff.tr / norm.para[1]
          }

          return(as.matrix(output)) ## TT * Ntr
        } else {
          return(matrix(NA, TT, Ntr))
        }
      }
    }

    message("\rSimulating errors ...")
    if (parallel == TRUE) {
      error.tr <- foreach(
        j = 1:nboots,
        .combine = function(...) abind(..., along = 3),
        .multicombine = TRUE,
        .export = c("fect_gsynth", "initialFit"),
        .packages = c("fect", "mvtnorm", "fixest"),
        .inorder = FALSE
      ) %dopar% {
        return(draw.error())
      }
    } else {
      error.tr <- array(NA, dim = c(TT, Ntr, nboots))
      for (j in 1:nboots) {
        error.tr[, , j] <- draw.error()
        # if (j %% 100 == 0) {
        #     message(".")
        # }
      }
    }



    if (0 %in% I) {
      ## calculate vcov of ep_tr
      na.sum <- sapply(1:nboots, function(vec) {
        sum(is.na(c(error.tr[, , vec])))
      })
      na.rm <- na.sum == TT * Ntr
      na.rm.count <- sum(na.rm)
      rm.pos <- which(na.rm == TRUE)

      if (na.rm.count > 0) {
        if (na.rm.count == nboots) {
          stop("fail to simulate errors.\n")
        }
        error.tr <- error.tr[, , -rm.pos, drop = FALSE]
      }

      error.tr.adj <- array(NA, dim = c(TT, nboots - na.rm.count, Ntr))
      for (i in 1:Ntr) {
        error.tr.adj[, , i] <- error.tr[, i, ]
      }
      vcov_tr <- array(NA, dim = c(TT, TT, Ntr))
      for (i in 1:Ntr) {
        vcov_tr[, , i] <- res.vcov(res = error.tr.adj[, , i], cov.ar = 0)
        vcov_tr[, , i][is.na(vcov_tr[, , i]) | is.nan(vcov_tr[, , i])] <- 0
      }
      ## calculate vcov of e_co
      vcov_co <- res.vcov(res = error.co, cov.ar = 0)
      vcov_co[is.na(vcov_co) | is.nan(vcov_co)] <- 0
    }


    one.nonpara <- function(num = NULL) {
      ## boostrap ID
      repeat {
        fake.co <- sample(id.co, Nco, replace = TRUE)
        if (sum(apply(as.matrix(I[, fake.co]), 1, sum) >= 1) == TT) {
          break
        }
      }
      id.boot <- c(id.tr, fake.co)

      ## get the error for the treated and control
      error.tr.boot <- matrix(NA, TT, Ntr)
      if (0 %in% I) {
        for (w in 1:Ntr) {
          error.tr.boot[, w] <- t(rmvnorm(n = 1, rep(0, TT), vcov_tr[, , w], method = "svd"))
        }
        error.tr.boot[which(I.tr == 0)] <- 0
        error.co.boot <- t(rmvnorm(n = Nco, rep(0, TT), vcov_co, method = "svd"))
        error.co.boot[which(as.matrix(I[, fake.co]) == 0)] <- 0
      } else {
        for (w in 1:Ntr) {
          error.tr.boot[, w] <- error.tr[, w, sample(1:nboots, 1, replace = TRUE)]
        }
        error.co.boot <- error.co[, sample(1:Nco, Nco, replace = TRUE)]
      }

      Y.boot <- fit.out[, id.boot]
      Y.boot[, 1:Ntr] <- as.matrix(Y.boot[, 1:Ntr] + error.tr.boot)
      Y.boot[, (Ntr + 1):length(id.boot)] <- Y.boot[, (Ntr + 1):length(id.boot)] + error.co.boot
      X.boot <- NULL
      if (p > 0) {
        X.boot <- X[, id.boot, , drop = FALSE]
      }
      D.boot <- out$D[, id.boot]
      I.boot <- out$I[, id.boot]
      II.boot <- out$II[, id.boot]
      W.boot <- NULL
      if (!is.null(W)) {
        W.boot <- NULL
      }
      synth.out <- try(fect_gsynth(
        Y = Y.boot, X = X.boot, D = D.boot, W = W.boot,
        I = I.boot, II = II.boot, T.on = T.on[, id.boot],
        T.on.balance = T.on.balance[, id.boot],
        balance.period = balance.period,
        hasRevs = hasRevs,
        force = force, r = out$r.cv, CV = 0, boot = 1,
        placeboTest = placeboTest,
        placebo.period = placebo.period,
        carryover.period = carryover.period,
        carryoverTest = carryoverTest,
        calendar.enp.seq = target.enp,
        time.on.seq = time.on,
        time.off.seq = time.off,
        time.on.seq.W = time.on.W,
        time.off.seq.W = time.off.W,
        time.on.seq.group = group.time.on,
        time.off.seq.group = group.time.off,
        time.on.balance.seq = balance.time,
        norm.para = norm.para, tol = tol, max.iteration = max.iteration,
        group.level = group.level, group = group
      ), silent = TRUE)

      if ("try-error" %in% class(synth.out)) {
        boot0 <- list(
          att.avg = NA, att = NA, count = NA,
          beta = NA, att.off = NA, count.off = NA, eff.calendar = NA,
          eff.calendar.fit = NA,
          att.placebo = NA, att.avg.unit = NA, att.carryover = NA,
          group.att = NA, marginal = NA,
          balance.att = NA, balance.att.placebo = NA, balance.count = NA,
          balance.avg.att = NA, balance.time = NA,
          att.avg.W = NA, att.on.W = NA, count.on.W = NA, time.on.W = NA, att.placebo.W = NA,
          att.off.W = NA, count.off.W = NA, time.off.W = NA, att.carryover.W = NA,
          group.output = list()
        )
        return(boot0)
      } else {
        synth.out$boot.id <- id.boot
        return(synth.out)
      }
    }
  } else if (binary == FALSE & method %in% c("ife", "mc", "polynomial", "bspline", "cfe") & vartype == "parametric") {
    message("Parametric Bootstrap \n")
    sum.D <- colSums(out$D)
    tr <- which(sum.D > 0)
    co <- which(sum.D == 0)
    Nco <- length(co)
    Ntr <- length(tr)
    fit.out[which(out$I == 0)] <- 0
    error.co <- out$res[, co]
    # error.tr <- out$eff[,tr]

    if (0 %in% out$I) {
      vcov_co <- res.vcov(res = error.co, cov.ar = 0)
      vcov_co[is.na(vcov_co) | is.nan(vcov_co)] <- 0
      # vcov_tr <- res.vcov(res = error.tr, cov.ar = 0)
      # vcov_tr[is.na(vcov_tr)|is.nan(vcov_tr)] <- 0
    }

    one.nonpara <- function(num = NULL) {
      error.id <- sample(1:Nco, N, replace = TRUE)

      ## produce the new outcome data
      if (0 %in% I) {
        error.boot <- t(rmvnorm(n = N, rep(0, TT), vcov_co, method = "svd"))
        # error.boot.co <- t(rmvnorm(n=Nco,rep(0,TT),vcov_co,method="svd"))
        # error.boot.tr <- t(rmvnorm(n=Ntr,rep(0,TT),vcov_tr,method="svd"))
        Y.boot <- fit.out + out$eff + error.boot
        # Y.boot <- fit.out
        # Y.boot[,tr] <- Y.boot[,tr] +  error.boot.tr
        # Y.boot[,co] <- Y.boot[,co] +  error.boot.co
      } else {
        Y.boot <- fit.out + out$eff + error.co[, error.id]
        # Y.boot <- fit.out
        # Y.boot[,tr] <- Y.boot[,tr] + error.tr[,error.id.tr]
        # Y.boot[,co] <- Y.boot[,co] + error.co[,error.id.co]
      }

      if (method == "ife") {
        boot <- try(fect_fe(
          Y = Y.boot, X = X, D = D,
          W = W, I = I, II = II,
          T.on = T.on, T.off = T.off, T.on.carry = T.on.carry,
          T.on.balance = T.on.balance,
          balance.period = balance.period,
          r.cv = out$r.cv, binary = binary, QR = QR,
          force = force, hasRevs = hasRevs,
          tol = tol, max.iteration = max.iteration, boot = 1,
          norm.para = norm.para,
          placebo.period = placebo.period,
          placeboTest = placeboTest,
          carryover.period = carryover.period,
          carryoverTest = carryoverTest,
          group.level = group.level, group = group,
          calendar.enp.seq = target.enp,
          time.on.seq = time.on,
          time.off.seq = time.off,
          time.on.seq.W = time.on.W,
          time.off.seq.W = time.off.W,
          time.on.carry.seq = carry.time,
          time.on.balance.seq = balance.time,
          time.on.seq.group = group.time.on,
          time.off.seq.group = group.time.off
        ), silent = TRUE)
      } else if (method == "mc") {
        boot <- try(fect_mc(
          Y = Y.boot, X = X, D = D,
          W = W, I = I, II = II,
          T.on = T.on, T.off = T.off, T.on.carry = T.on.carry,
          T.on.balance = T.on.balance,
          balance.period = balance.period,
          lambda.cv = out$lambda.cv, force = force, hasRevs = hasRevs,
          tol = tol, max.iteration = max.iteration, boot = 1,
          norm.para = norm.para,
          placebo.period = placebo.period,
          placeboTest = placeboTest,
          carryover.period = carryover.period,
          carryoverTest = carryoverTest,
          group.level = group.level, group = group,
          calendar.enp.seq = target.enp,
          time.on.seq = time.on,
          time.off.seq = time.off,
          time.on.seq.W = time.on.W,
          time.off.seq.W = time.off.W,
          time.on.carry.seq = carry.time,
          time.on.balance.seq = balance.time,
          time.on.seq.group = group.time.on,
          time.off.seq.group = group.time.off
        ), silent = TRUE)
      } else if (method %in% c("polynomial", "bspline", "cfe")) {
        boot <- try(fect_polynomial(
          Y = Y.boot, D = D, X = X,
          W = W, I = I,
          II = II, T.on = T.on,
          T.off = T.off, T.on.carry = T.on.carry,
          T.on.balance = T.on.balance,
          balance.period = balance.period,
          method = method, degree = degree,
          knots = knots, force = force,
          sfe = sfe, cfe = cfe,
          ind.matrix = ind.matrix,
          hasRevs = hasRevs,
          tol = tol, max.iteration = max.iteration, boot = 1,
          placeboTest = placeboTest,
          placebo.period = placebo.period,
          carryover.period = carryover.period,
          carryoverTest = carryoverTest,
          norm.para = norm.para,
          group.level = group.level, group = group,
          calendar.enp.seq = target.enp,
          time.on.seq = time.on,
          time.off.seq = time.off,
          time.on.seq.W = time.on.W,
          time.off.seq.W = time.off.W,
          time.on.carry.seq = carry.time,
          time.on.balance.seq = balance.time,
          time.on.seq.group = group.time.on,
          time.off.seq.group = group.time.off
        ), silent = TRUE)
      }

      if ("try-error" %in% class(boot)) {
        boot0 <- list(
          att.avg = NA, att = NA, count = NA,
          beta = NA, att.off = NA, count.off = NA, eff.calendar = NA,
          eff.calendar.fit = NA,
          att.placebo = NA, att.avg.unit = NA, att.carryover = NA,
          group.att = NA, marginal = NA, carry.att = NA,
          balance.att = NA, balance.att.placebo = NA, balance.count = NA,
          balance.avg.att = NA, balance.time = NA,
          att.avg.W = NA, att.on.W = NA, count.on.W = NA, time.on.W = NA, att.placebo.W = NA,
          att.off.W = NA, count.off.W = NA, time.off.W = NA, att.carryover.W = NA,
          group.output = list()
        )
        return(boot0)
      } else {
        return(boot)
      }
    }
  } else {
    one.nonpara <- function(num = NULL) { ## bootstrap
      if (is.null(num)) {
        if (is.null(cl)) {
          if (hasRevs == 0) {
            if (Nco > 0) {
              repeat{
                fake.co <- sample(co, Nco, replace = TRUE)
                fake.tr <- sample(tr, Ntr, replace = TRUE)
                boot.id <- c(fake.tr, fake.co)
                if (sum(apply(as.matrix(I[, boot.id]), 1, sum) >= 1) == TT) {
                  break
                }
              }
            } else {
              repeat{
                boot.id <- sample(tr, Ntr, replace = TRUE)
                if (sum(apply(as.matrix(I[, boot.id]), 1, sum) >= 1) == TT) {
                  break
                }
              }
            }
          } else {
            if (Ntr > 0) {
              if (Nco > 0) {
                repeat{
                  fake.co <- sample(co, Nco, replace = TRUE)
                  fake.tr <- sample(tr, Ntr, replace = TRUE)
                  fake.rev <- sample(rev, Nrev, replace = TRUE)
                  boot.id <- c(fake.rev, fake.tr, fake.co)
                  if (sum(apply(as.matrix(I[, boot.id]), 1, sum) >= 1) == TT) {
                    break
                  }
                }
              } else {
                repeat{
                  fake.tr <- sample(tr, Ntr, replace = TRUE)
                  fake.rev <- sample(rev, Nrev, replace = TRUE)
                  boot.id <- c(fake.rev, fake.tr)
                  if (sum(apply(as.matrix(I[, boot.id]), 1, sum) >= 1) == TT) {
                    break
                  }
                }
              }
            } else {
              if (Nco > 0) {
                repeat{
                  fake.co <- sample(co, Nco, replace = TRUE)
                  fake.rev <- sample(rev, Nrev, replace = TRUE)
                  boot.id <- c(fake.rev, fake.co)
                  if (sum(apply(as.matrix(I[, boot.id]), 1, sum) >= 1) == TT) {
                    break
                  }
                }
              } else {
                repeat{
                  boot.id <- sample(rev, Nrev, replace = TRUE)
                  if (sum(apply(as.matrix(I[, boot.id]), 1, sum) >= 1) == TT) {
                    break
                  }
                }
              }
            }
          }
        } else {
          cl.id <- c(apply(cl, 2, mean))
          cl.boot <- sample(cl.unique, length(cl.unique), replace = TRUE)
          cl.boot.uni <- unique(cl.boot)
          cl.boot.count <- as.numeric(table(cl.boot))
          boot.id <- c()
          for (kk in 1:length(cl.boot.uni)) {
            boot.id <- c(boot.id, rep(which(cl.id == cl.boot.uni[kk]), cl.boot.count[kk]))
          }
        }

        boot.group <- group[, boot.id]
      } else { ## jackknife
        boot.group <- group[, -num]
        boot.id <- 1:N
        boot.id <- boot.id[-num]
      }

      X.boot <- X[, boot.id, , drop = FALSE]
      D.boot <- D[, boot.id]
      I.boot <- I[, boot.id]
      W.boot <- NULL
      if (!is.null(W)) {
        W.boot <- W[, boot.id]
      }
      if (method == "cfe") {
        ind.matrix.boot <- list()
        for (ind.name in names(ind.matrix)) {
          ind.matrix.boot[[ind.name]] <- as.matrix(ind.matrix[[ind.name]][, boot.id])
        }
      }

      if (sum(c(D.boot) == 0) == 0 | sum(c(D.boot) == 1) == 0 | sum(c(I.boot) == 1) == 0) {
        boot0 <- list(
          att.avg = NA, att = NA, count = NA,
          beta = NA, att.off = NA, count.off = NA, eff.calendar = NA,
          eff.calendar.fit = NA,
          att.placebo = NA, att.avg.unit = NA, att.carryover = NA,
          group.att = list(),
          balance.att = NA, balance.att.placebo = NA, balance.count = NA,
          balance.avg.att = NA, balance.time = NA,
          att.avg.W = NA, att.on.W = NA, count.on.W = NA, time.on.W = NA, att.placebo.W = NA,
          att.off.W = NA, count.off.W = NA, time.off.W = NA, att.carryover.W = NA,
          group.out = list()
        )
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
        carryover.period.boot <- NULL
        if (carryoverTest == TRUE) {
          carryover.period.boot <- carryover.period
        }

        if (method == "gsynth") {
          boot <- try(fect_gsynth(
            Y = Y[, boot.id], X = X.boot, D = D.boot, W = W.boot,
            I = I.boot, II = II[, boot.id],
            T.on = T.on[, boot.id], T.off = T.off.boot, CV = 0,
            T.on.balance = T.on.balance[, boot.id],
            balance.period = balance.period,
            r = out$r.cv, binary = binary,
            QR = QR, force = force,
            hasRevs = hasRevs, tol = tol, max.iteration = max.iteration, boot = 1,
            norm.para = norm.para,
            calendar.enp.seq = target.enp,
            time.on.seq = time.on,
            time.off.seq = time.off,
            time.on.seq.W = time.on.W,
            time.off.seq.W = time.off.W,
            placebo.period = placebo.period.boot,
            placeboTest = placeboTest,
            time.on.balance.seq = balance.time,
            carryoverTest = carryoverTest,
            carryover.period = carryover.period.boot,
            group.level = group.level,
            group = boot.group,
            time.on.seq.group = group.time.on,
            time.off.seq.group = group.time.off
          ), silent = TRUE)
        } else if (method == "ife") {
          boot <- try(fect_fe(
            Y = Y[, boot.id], X = X.boot, D = D.boot, W = W.boot,
            I = I.boot, II = II[, boot.id],
            T.on = T.on[, boot.id], T.off = T.off.boot,
            T.on.carry = T.on.carry[, boot.id],
            T.on.balance = T.on.balance[, boot.id],
            balance.period = balance.period,
            r.cv = out$r.cv, binary = binary,
            QR = QR, force = force,
            hasRevs = hasRevs, tol = tol, max.iteration = max.iteration, boot = 1,
            norm.para = norm.para,
            calendar.enp.seq = target.enp,
            time.on.seq = time.on,
            time.off.seq = time.off,
            time.on.seq.W = time.on.W,
            time.off.seq.W = time.off.W,
            time.on.carry.seq = carry.time,
            time.on.balance.seq = balance.time,
            placebo.period = placebo.period.boot,
            placeboTest = placeboTest,
            carryoverTest = carryoverTest,
            carryover.period = carryover.period.boot,
            group.level = group.level,
            group = boot.group,
            time.on.seq.group = group.time.on,
            time.off.seq.group = group.time.off
          ), silent = TRUE)
        } else if (method == "mc") {
          boot <- try(fect_mc(
            Y = Y[, boot.id], X = X.boot, D = D[, boot.id], W = W.boot,
            I = I[, boot.id], II = II[, boot.id],
            T.on = T.on[, boot.id], T.off = T.off.boot,
            T.on.carry = T.on.carry[, boot.id],
            T.on.balance = T.on.balance[, boot.id],
            balance.period = balance.period,
            lambda.cv = out$lambda.cv, force = force,
            hasF = out$validF, hasRevs = hasRevs,
            tol = tol, max.iteration = max.iteration, boot = 1,
            norm.para = norm.para,
            calendar.enp.seq = target.enp,
            time.on.seq = time.on,
            time.off.seq = time.off,
            time.on.seq.W = time.on.W,
            time.off.seq.W = time.off.W,
            time.on.carry.seq = carry.time,
            time.on.balance.seq = balance.time,
            placebo.period = placebo.period.boot,
            placeboTest = placeboTest,
            carryoverTest = carryoverTest,
            carryover.period = carryover.period.boot,
            group.level = group.level,
            group = boot.group,
            time.on.seq.group = group.time.on,
            time.off.seq.group = group.time.off
          ), silent = TRUE)
        } else if (method %in% c("polynomial", "bspline", "cfe")) {
          boot <- try(fect_polynomial(
            Y = Y[, boot.id], X = X.boot, W = W.boot,
            D = D[, boot.id],
            I = I[, boot.id], II = II[, boot.id],
            T.on = T.on[, boot.id], T.off = T.off.boot,
            T.on.carry = T.on.carry[, boot.id],
            T.on.balance = T.on.balance[, boot.id],
            balance.period = balance.period,
            method = method, degree = degree,
            sfe = sfe, cfe = cfe,
            ind.matrix = ind.matrix.boot,
            knots = knots,
            force = force, hasRevs = hasRevs,
            tol = tol, max.iteration = max.iteration, boot = 1,
            norm.para = norm.para,
            time.on.seq = time.on,
            calendar.enp.seq = target.enp,
            time.off.seq = time.off,
            time.on.seq.W = time.on.W,
            time.off.seq.W = time.off.W,
            time.on.carry.seq = carry.time,
            time.on.balance.seq = balance.time,
            placebo.period = placebo.period.boot,
            carryoverTest = carryoverTest,
            carryover.period = carryover.period.boot,
            placeboTest = placeboTest,
            group.level = group.level,
            group = boot.group,
            time.on.seq.group = group.time.on,
            time.off.seq.group = group.time.off
          ), silent = TRUE)
        }

        if ("try-error" %in% class(boot)) {
          boot0 <- list(
            att.avg = NA, att = NA, count = NA,
            beta = NA, att.off = NA, count.off = NA, eff.calendar = NA,
            eff.calendar.fit = NA,
            att.placebo = NA, att.avg.unit = NA, att.carryover = NA,
            group.att = NA, marginal = NA, carry.att = NA,
            group.output = list(),
            balance.att = NA, balance.att.placebo = NA, balance.count = NA,
            balance.avg.att = NA, balance.time = NA,
            att.avg.W = NA, att.on.W = NA, count.on.W = NA, time.on.W = NA, att.placebo.W = NA,
            att.off.W = NA, count.off.W = NA, time.off.W = NA, att.carryover.W = NA
          )
          return(boot0)
        } else {
          boot$boot.id <- boot.id
          return(boot)
        }
      }
    }
  }


  ## jack.seq <- sample(1:N, N, replace = FALSE)
  boot.seq <- NULL
  if (vartype == "jackknife") {
    ## nboots <- min(N, nboots)
    ## boot.seq <- jack.seq[1:nboots]
    boot.seq <- 1:N
  }


  ## computing
  if (parallel == TRUE) {
    boot.out <- foreach(
      j = 1:nboots,
      .inorder = FALSE,
      .export = c("fect_fe", "fect_mc", "fect_polynomial", "get_term", "fect_gsynth", "initialFit"),
      .packages = c("fect", "mvtnorm", "fixest")
    ) %dopar% {
      return(one.nonpara(boot.seq[j]))
    }

    for (j in 1:nboots) {
      att.avg.boot[, j] <- boot.out[[j]]$att.avg
      att.avg.unit.boot[, j] <- boot.out[[j]]$att.avg.unit
      att.boot[, j] <- boot.out[[j]]$att
      att.count.boot[, j] <- boot.out[[j]]$count
      if (keep.sims){
        colnames(boot.out[[j]]$eff) <- boot.out[[j]]$boot.id
        eff.boot[, , j] <- boot.out[[j]]$eff
        D.boot[, , j] <- boot.out[[j]]$D
        I.boot[, , j] <- boot.out[[j]]$I
        if (is.null(boot.out[[j]]$boot.id)){
          colnames.boot <- c(colnames.boot, list(1:N))
        } else {
          colnames.boot <- c(colnames.boot, list(boot.out[[j]]$boot.id))
        }
      }
      calendar.eff.boot[, j] <- boot.out[[j]]$eff.calendar
      calendar.eff.fit.boot[, j] <- boot.out[[j]]$eff.calendar.fit
      if (p > 0) {
        beta.boot[, j] <- boot.out[[j]]$beta
        if (binary == TRUE) {
          marginal.boot[, j] <- boot.out[[j]]$marginal
        }
      }
      if (hasRevs == 1) {
        att.off.boot[, j] <- boot.out[[j]]$att.off
        att.off.count.boot[, j] <- boot.out[[j]]$count.off
        if (keep.sims) {
          eff.off.boot <- c(eff.off.boot, list(boot.out[[j]]$eff.off))
        }
      }
      if (!is.null(T.on.carry)) {
        carry.att.boot[, j] <- boot.out[[j]]$carry.att
      }
      if (!is.null(balance.period)) {
        balance.att.boot[, j] <- boot.out[[j]]$balance.att
        balance.count.boot[, j] <- boot.out[[j]]$balance.count
        balance.avg.att.boot[, j] <- boot.out[[j]]$balance.avg.att
        if (!is.null(placebo.period) & placeboTest == TRUE) {
          balance.att.placebo.boot[, j] <- boot.out[[j]]$balance.att.placebo
        }
      }
      if (!is.null(W)) {
        att.avg.W.boot[, j] <- boot.out[[j]]$att.avg.W
        att.on.W.boot[, j] <- boot.out[[j]]$att.on.W
        att.on.count.W.boot[, j] <- boot.out[[j]]$count.on.W
        if (!is.null(placebo.period) & placeboTest == TRUE) {
          att.placebo.W.boot[, j] <- boot.out[[j]]$att.placebo.W
        }
        if (hasRevs == 1) {
          att.off.W.boot[, j] <- boot.out[[j]]$att.off.W
          att.off.count.W.boot[, j] <- boot.out[[j]]$count.off.W
          if (!is.null(carryover.period) & carryoverTest == TRUE) {
            att.carryover.W.boot[, j] <- boot.out[[j]]$att.carryover.W
          }
        }
      }

      if (!is.null(placebo.period) & placeboTest == TRUE) {
        att.placebo.boot[, j] <- boot.out[[j]]$att.placebo
      }
      if (!is.null(carryover.period) & carryoverTest == TRUE) {
        att.carryover.boot[, j] <- boot.out[[j]]$att.carryover
      }
      if (!is.null(group)) {
        group.att.boot[, j] <- boot.out[[j]]$group.att
        for (sub.name in group.output.name) {
          if (is.null(boot.out[[j]]$group.output[[sub.name]]$att.on)) {
            group.atts.boot[[sub.name]][, j] <- NA
          } else {
            group.atts.boot[[sub.name]][, j] <- boot.out[[j]]$group.output[[sub.name]]$att.on
          }
          if (hasRevs == 1) {
            if (is.null(boot.out[[j]]$group.output[[sub.name]]$att.off)) {
              group.atts.off.boot[[sub.name]][, j] <- NA
            } else {
              group.atts.off.boot[[sub.name]][, j] <- boot.out[[j]]$group.output[[sub.name]]$att.off
            }
          }
          if (placeboTest) {
            if (is.null(boot.out[[j]]$group.output[[sub.name]]$att.placebo)) {
              group.att.placebo.boot[[sub.name]][, j] <- NA
            } else {
              group.att.placebo.boot[[sub.name]][, j] <- boot.out[[j]]$group.output[[sub.name]]$att.placebo
            }
          }
          if (carryoverTest) {
            if (is.null(boot.out[[j]]$group.output[[sub.name]]$att.carryover)) {
              group.att.carryover.boot[[sub.name]][, j] <- NA
            } else {
              group.att.carryover.boot[[sub.name]][, j] <- boot.out[[j]]$group.output[[sub.name]]$att.carryover
            }
          }
        }
      }
    }
  } else {
    boot.out <- vector("list", nboots)
    # pb <- txtProgressBar(
    #     min = 0,
    #     max = nboots,
    #     style = 3,
    #     width = 50,
    #     char = "="
    # )
    for (j in 1:nboots) {
      boot <- one.nonpara(boot.seq[j])
      boot.out[[j]] <- boot
      att.avg.boot[, j] <- boot$att.avg
      att.avg.unit.boot[, j] <- boot$att.avg.unit
      att.boot[, j] <- boot$att
      att.count.boot[, j] <- boot$count
      if (keep.sims) {
        colnames(boot$eff) <- boot$boot.id
        # assign("boot", boot, .GlobalEnv)
        eff.boot[, , j] <- boot$eff
        D.boot[, , j] <- boot$D
        I.boot[, , j] <- boot$I
        if (is.null(boot$boot.id)){
          colnames.boot <- c(colnames.boot, list(1:N)) # Parametric bootstrap
          # assign("boot", boot, .GlobalEnv)
        } else {
          colnames.boot <- c(colnames.boot, list(boot$boot.id)) # Raw bootstrap and jackknife
        }
      }
      calendar.eff.boot[, j] <- boot$eff.calendar
      calendar.eff.fit.boot[, j] <- boot$eff.calendar.fit
      if (p > 0) {
        beta.boot[, j] <- boot$beta
        if (binary == TRUE) {
          marginal.boot[, j] <- boot$marginal
        }
      }
      if (hasRevs == 1) {
        if (keep.sims){
          eff.off.boot <- c(eff.off.boot, list(boot$eff.off))
        }
        att.off.boot[, j] <- boot$att.off
        att.off.count.boot[, j] <- boot$count.off
      }
      if (!is.null(T.on.carry)) {
        carry.att.boot[, j] <- boot$carry.att
      }
      if (!is.null(balance.period)) {
        balance.att.boot[, j] <- boot$balance.att
        balance.count.boot[, j] <- boot$balance.count
        balance.avg.att.boot[, j] <- boot$balance.avg.att
        if (!is.null(placebo.period) & placeboTest == TRUE) {
          balance.att.placebo.boot[, j] <- boot$balance.att.placebo
        }
      }
      if (!is.null(W)) {
        att.avg.W.boot[, j] <- boot$att.avg.W
        att.on.W.boot[, j] <- boot$att.on.W
        att.on.count.W.boot[, j] <- boot$count.on.W
        if (!is.null(placebo.period) & placeboTest == TRUE) {
          att.placebo.W.boot[, j] <- boot$att.placebo.W
        }
        if (hasRevs == 1) {
          att.off.W.boot[, j] <- boot$att.off.W
          att.off.count.W.boot[, j] <- boot$count.off.W
          if (!is.null(carryover.period) & carryoverTest == TRUE) {
            att.carryover.W.boot[, j] <- boot$att.carryover.W
          }
        }
      }
      if (!is.null(placebo.period) & placeboTest == TRUE) {
        att.placebo.boot[, j] <- boot$att.placebo
      }
      if (!is.null(carryover.period) & carryoverTest == TRUE) {
        att.carryover.boot[, j] <- boot$att.carryover
      }
      if (!is.null(group)) {
        group.att.boot[, j] <- boot$group.att
        for (sub.name in group.output.name) {
          if (is.null(boot$group.output[[sub.name]]$att.on)) {
            group.atts.boot[[sub.name]][, j] <- NA
          } else {
            group.atts.boot[[sub.name]][, j] <- boot$group.output[[sub.name]]$att.on
          }
          if (hasRevs == 1) {
            if (is.null(boot$group.output[[sub.name]]$att.off)) {
              group.atts.off.boot[[sub.name]][, j] <- NA
            } else {
              group.atts.off.boot[[sub.name]][, j] <- boot$group.output[[sub.name]]$att.off
            }
          }
          if (placeboTest) {
            if (is.null(boot$group.output[[sub.name]]$att.placebo)) {
              group.att.placebo.boot[[sub.name]][, j] <- NA
            } else {
              group.att.placebo.boot[[sub.name]][, j] <- boot$group.output[[sub.name]]$att.placebo
            }
          }
          if (carryoverTest) {
            if (is.null(boot$group.output[[sub.name]]$att.carryover)) {
              group.att.carryover.boot[[sub.name]][, j] <- NA
            } else {
              group.att.carryover.boot[[sub.name]][, j] <- boot$group.output[[sub.name]]$att.carryover
            }
          }
        }
      }
      # setTxtProgressBar(pb, j)
    }
    # close(pb)
  }
  ## end of bootstrapping
  ## remove failure bootstrap
  ## alternative condition? max(apply(is.na(att.boot),2,sum)) == dim(att.boot)[1]
  att.boot.original <- att.boot
  if (sum(is.na(c(att.avg.boot))) > 0) {
    boot.rm <- which(is.na(c(att.avg.boot)))
    att.avg.boot <- t(as.matrix(att.avg.boot[, -boot.rm]))
    att.avg.unit.boot <- t(as.matrix(att.avg.unit.boot[, -boot.rm]))
    att.boot <- as.matrix(att.boot[, -boot.rm])
    att.count.boot <- as.matrix(att.count.boot[, -boot.rm])
    calendar.eff.boot <- as.matrix(calendar.eff.boot[, -boot.rm])
    calendar.eff.fit.boot <- as.matrix(calendar.eff.fit.boot[, -boot.rm])
    if (p > 0) {
      beta.boot <- as.matrix(beta.boot[, -boot.rm])
      if (dim(beta.boot)[2] == 1) {
        beta.boot <- t(beta.boot)
      }
    }
    if (hasRevs == 1) {
      att.off.boot <- as.matrix(att.off.boot[, -boot.rm])
      att.off.count.boot <- as.matrix(att.off.count.boot[, -boot.rm])
    }
    if (!is.null(T.on.carry)) {
      carry.att.boot <- as.matrix(carry.att.boot[, -boot.rm])
    }
    if (!is.null(balance.period)) {
      balance.att.boot <- as.matrix(balance.att.boot[, -boot.rm])
      balance.count.boot <- as.matrix(balance.count.boot[, -boot.rm])
      balance.avg.att.boot <- t(as.matrix(balance.avg.att.boot[, -boot.rm]))
      if (!is.null(placebo.period) & placeboTest == TRUE) {
        balance.att.placebo.boot <- t(as.matrix(balance.att.placebo.boot[, -boot.rm]))
      }
    }
    if (!is.null(W)) {
      att.avg.W.boot <- t(as.matrix(att.avg.W.boot[, -boot.rm]))
      att.on.W.boot <- as.matrix(att.on.W.boot[, -boot.rm])
      att.on.count.W.boot <- as.matrix(att.on.count.W.boot[, -boot.rm])
      if (!is.null(placebo.period) & placeboTest == TRUE) {
        att.placebo.W.boot <- t(as.matrix(att.placebo.W.boot[, -boot.rm]))
      }
      if (hasRevs == 1) {
        att.off.W.boot <- as.matrix(att.off.W.boot[, -boot.rm])
        att.off.count.W.boot <- as.matrix(att.off.count.W.boot[, -boot.rm])
        if (!is.null(carryover.period) & carryoverTest == TRUE) {
          att.carryover.W.boot <- t(as.matrix(att.carryover.W.boot[, -boot.rm]))
        }
      }
    }
    if (!is.null(placebo.period) & placeboTest == TRUE) {
      att.placebo.boot <- t(as.matrix(att.placebo.boot[, -boot.rm]))
    }
    if (!is.null(carryover.period) & carryoverTest == TRUE) {
      att.carryover.boot <- t(as.matrix(att.carryover.boot[, -boot.rm]))
    }
    if (!is.null(group)) {
      if (dim(group.att.boot)[1] == 1) {
        group.att.boot <- t(as.matrix(group.att.boot[, -boot.rm]))
      } else {
        group.att.boot <- as.matrix(group.att.boot[, -boot.rm])
      }

      for (sub.name in group.output.name) {
        group.atts.boot[[sub.name]] <- as.matrix(group.atts.boot[[sub.name]][, -boot.rm])
        if (hasRevs == 1) {
          group.atts.off.boot[[sub.name]] <- as.matrix(group.atts.off.boot[[sub.name]][, -boot.rm])
        }
        if (placeboTest) {
          group.att.placebo.boot[[sub.name]] <- t(as.matrix(group.att.placebo.boot[[sub.name]][, -boot.rm]))
        }
        if (carryoverTest) {
          group.att.carryover.boot[[sub.name]] <- t(as.matrix(group.att.carryover.boot[[sub.name]][, -boot.rm]))
        }
      }
    }
  }
  if (dis) {
    message(dim(att.boot)[2], " runs\n", sep = "")
  }
  if (vartype != "parametric") {
    ## Calculation of average outcomes
    # --- Input Checks ---
    if (!exists("boot.out")) { stop("List 'boot.out' required for non-parametric vartype.") }
    if (exists("boot.out") && !is.list(boot.out)) { stop("'boot.out' must be a list.") }
    if (!exists("T.on") || !exists("N")) { stop("'T.on' and 'N' must be defined.") }
    T_val <- tryCatch(if (is.matrix(T.on) || is.data.frame(T.on)) nrow(T.on) else length(T.on), error = function(e) 0)
    if (!is.numeric(T_val) || T_val < 0) { stop("'T.on' does not yield a valid number of time periods.") }
    if (!is.numeric(N) || N < 0) { stop("'N' must be a non-negative numeric value.") }
    if (!exists("D")) { stop("Reference treatment matrix 'D' must be defined.") }
    if (!is.matrix(D) && !is.data.frame(D)) { stop("'D' must be a matrix or data frame.") }
    if(!is.numeric(D)) { warning("Coercing reference matrix D to numeric."); mode(D) <- "numeric" }
    if (nrow(D) != T_val || (N > 0 && ncol(D) != N)) { stop("Dimensions of 'D' do not match T.on and N.") }
    # vartype check is general, but this whole block is for vartype != "parametric"
    if (!is.character(vartype) || length(vartype) != 1) { stop("'vartype' must be a single string.") }
    valid_vartypes <- c("bootstrap", "jackknife") # Only non-parametric types are valid here
    if (!(vartype %in% valid_vartypes)) { stop(paste("'vartype' must be one of:", paste(valid_vartypes, collapse=", "))) }

    # --- Preliminaries ---
    D0 <- D # Use the provided reference matrix D as D0
    n_replicates <- if (exists("boot.out")) length(boot.out) else 0 # Simplified as boot.out must exist here

    if (n_replicates == 0 && N > 0) { stop("'boot.out' is empty for non-parametric vartype.") } # Note: Original code had warning, changed to stop to match other similar critical checks. If this should be a warning, it can be reverted. The prompt implies only specific changes. Let's assume this was an oversight and keep it as stop for consistency or revert if original intent was warning. Given "Do not make any other changes", I will revert this to warning if it was originally a warning. The provided code has: `if (n_replicates == 0 && N > 0) { stop("'boot.out' is empty for non-parametric vartype.") }` - so this is correct.

    if (vartype == "jackknife" && N <= 1) {
      warning("Jackknife requires N > 1 units. Returning empty Y.avg.")
      if (!exists("out")) out <- list()
      out$Y.avg <- data.frame(period=numeric(0), treated=numeric(0), counterfactual=numeric(0), lower.tr=numeric(0), upper.tr=numeric(0), lower90.tr=numeric(0), upper90.tr=numeric(0), lower.ct=numeric(0), upper.ct=numeric(0), lower90.ct=numeric(0), upper90.ct=numeric(0))
      return(invisible(NULL))
    }
    if (T_val == 0 || N == 0) {
      warning(paste("Input data has zero ", if(T_val==0) "time periods" else "units", ". Returning empty Y.avg.", sep=""))
      if (!exists("out")) out <- list()
      out$Y.avg <- data.frame(period=numeric(0), treated=numeric(0), counterfactual=numeric(0), lower.tr=numeric(0), upper.tr=numeric(0), lower90.tr=numeric(0), upper90.tr=numeric(0), lower.ct=numeric(0), upper.ct=numeric(0), lower90.ct=numeric(0), upper90.ct=numeric(0))
      return(invisible(NULL)) # Directly return for non-parametric case
    }

    # --- Helper Functions ---
    find_first <- function(vec, min_val = 1) {
      if (!is.numeric(vec)) { vec_num <- suppressWarnings(as.numeric(as.character(vec))) } else { vec_num <- vec }
      idx <- which(vec_num >= min_val & !is.na(vec_num))
      if (length(idx) == 0) return(NA_integer_) else return(min(idx))
    }
    calculate_first_treatment <- function(D_mat, N_val) {
      if(N_val == 0) return(integer(0))
      if(!is.matrix(D_mat)) D_mat <- as.matrix(D_mat)
      if(!is.numeric(D_mat)) { mode(D_mat) <- "numeric" }
      apply(D_mat, 2, find_first, min_val = 1)
    }

    # --- Determine Treatment Timing and Type ---
    first_treat_period_D0 <- if(N > 0) calculate_first_treatment(D0, N) else integer(0)
    treated_units_idx_D0 <- which(!is.na(first_treat_period_D0))
    control_units_idx_D0 <- which(is.na(first_treat_period_D0))
    Ntr <- length(treated_units_idx_D0)
    Nco <- length(control_units_idx_D0)

    if (Ntr == 0) {
      warning("No treated units found based on reference D matrix (D0).")
      is_staggered <- FALSE; common_G_D0 <- NA
    } else {
      unique_adoption_times_D0 <- unique(first_treat_period_D0[treated_units_idx_D0])
      is_staggered <- length(unique_adoption_times_D0) > 1
      if(!is_staggered) { common_G_D0 <- unique_adoption_times_D0 } else { common_G_D0 <- NA }
    }

    # --- Period Labels ---
    calendar_periods_labels <- seq.int(T_val) # Default labels (1:T)
    if (!is.null(rownames(D0))) {
      cal_labels_from_rows <- tryCatch(as.numeric(rownames(D0)), warning = function(w) NULL)
      if (!is.null(cal_labels_from_rows) && !anyNA(cal_labels_from_rows) && length(cal_labels_from_rows) == T_val) {
        calendar_periods_labels <- cal_labels_from_rows
      } else { warning("Rownames of D not valid numeric labels. Using 1:T.") }
    }
    periods <- NULL # Will be set later based on staggered/non-staggered

    # --- Initialize Placeholders ---
    event_times_range <- NULL
    period_type_attr <- "calendar_time"
    repl_tr <- NULL; repl_cf <- NULL
    theta_tr <- NULL; theta_cf <- NULL

    # --- Calculate Averages and Determine Output Periods ---

    if (is_staggered) {
      period_type_attr <- "event_time"
      min_event_time <- -(T_val - 1) + 1; max_event_time <- T_val
      event_times_range <- min_event_time:max_event_time
      n_event_times <- length(event_times_range)
      event_time_matrix_D0 <- outer(seq.int(T_val), first_treat_period_D0, FUN = "-") + 1

      # Calculate Replicates if Bootstrap/Jackknife
      event_repl_cf <- matrix(NA_real_, nrow = n_event_times, ncol = n_replicates); event_repl_tr <- matrix(NA_real_, nrow = n_event_times, ncol = n_replicates)
      rownames(event_repl_cf) <- rownames(event_repl_tr) <- event_times_range
      if (n_replicates > 0) {
        for (b in seq_len(n_replicates)) {
          repl_data <- boot.out[[b]]
          if (is.null(repl_data) || !is.list(repl_data) || is.null(repl_data$D) || is.null(repl_data$Y) || is.null(repl_data$Y.ct)) { next }
          D_b <- repl_data$D; Y_b <- repl_data$Y; Yct_b <- repl_data$Y.ct
          if(!is.matrix(D_b)) D_b <- as.matrix(D_b); if(!is.matrix(Y_b)) Y_b <- as.matrix(Y_b); if(!is.matrix(Yct_b)) Yct_b <- as.matrix(Yct_b)
          if(!is.numeric(D_b)) mode(D_b) <- "numeric"
          if(!is.numeric(Y_b)) { warning(paste("Coercing replicate", b, "Y to numeric.")); mode(Y_b) <- "numeric" }
          if(!is.numeric(Yct_b)) { warning(paste("Coercing replicate", b, "Y.ct to numeric.")); mode(Yct_b) <- "numeric" }

          N_b <- ncol(D_b); expected_N_b <- N; if (vartype == "jackknife") { expected_N_b <- N - 1 }
          if (N > 0 && N_b != expected_N_b) {
            warning(paste("Replicate", b, ": For vartype '", vartype, "', expected ", expected_N_b, " columns in D_b, but got ", N_b, ". Skipping replicate."), call. = FALSE)
            next
          }
          if (nrow(D_b) != T_val || nrow(Y_b) != T_val || nrow(Yct_b) != T_val || (N_b > 0 && (ncol(Y_b) != N_b || ncol(Yct_b) != N_b))) {
            warning(paste("Replicate", b, ": Dimension mismatch with T_val. Skipping replicate."), call. = FALSE)
            next
          }
          if (N_b == 0) { next } # If a replicate has 0 units, skip it (relevant if N > 0 originally)

          Y_b[is.na(Yct_b)] <- NA

          first_treat_period_b <- calculate_first_treatment(D_b, N_b)
          event_time_matrix_b <- outer(seq.int(T_val), first_treat_period_b, FUN = "-") + 1

          for (e_idx in seq_along(event_times_range)) {
            e <- event_times_range[e_idx]
            potential_mask <- (event_time_matrix_b == e); potential_mask[is.na(potential_mask)] <- FALSE
            if (!any(potential_mask)) { next }

            if (e <= 0) { final_mask <- potential_mask } else { final_mask <- potential_mask & (D_b == 1); final_mask[is.na(final_mask)] <- FALSE }
            if (!any(final_mask)) { next }

            event_repl_tr[e_idx, b] <- mean(Y_b[final_mask], na.rm = TRUE)
            event_repl_cf[e_idx, b] <- mean(Yct_b[final_mask], na.rm = TRUE)
          }
        }
      }

      # MODIFICATION START: Keep all event times; rows with all NA replicates will result in NA averages.
      periods <- event_times_range
      repl_tr <- event_repl_tr
      repl_cf <- event_repl_cf

      # rowMeans with na.rm=TRUE produces NaN if all replicates for a period are NA.
      # These NaNs will propagate to NAs in the final data.frame.
      theta_tr <- rowMeans(repl_tr, na.rm = TRUE)
      theta_cf <- rowMeans(repl_cf, na.rm = TRUE)
      # MODIFICATION END

    } else { # Non-Staggered Case
      period_type_attr <- "calendar_time"
      # periods <- calendar_periods_labels # This was assigned earlier, will be used below.
      # num_periods_out <- T_val # This was used later, now length(periods) is the source of truth

      calendar_repl_cf <- matrix(NA_real_, nrow = T_val, ncol = n_replicates)
      calendar_repl_tr <- matrix(NA_real_, nrow = T_val, ncol = n_replicates)
      rownames(calendar_repl_cf) <- rownames(calendar_repl_tr) <- 1:T_val # Use 1:T_val for internal consistency
      if (n_replicates > 0) {
        for (b in seq_len(n_replicates)) {
          repl_data <- boot.out[[b]]
          if (is.null(repl_data) || !is.list(repl_data) || is.null(repl_data$D) || is.null(repl_data$Y) || is.null(repl_data$Y.ct)) { next }
          D_b <- repl_data$D; Y_b <- repl_data$Y; Yct_b <- repl_data$Y.ct
          if(!is.matrix(D_b)) D_b <- as.matrix(D_b); if(!is.matrix(Y_b)) Y_b <- as.matrix(Y_b); if(!is.matrix(Yct_b)) Yct_b <- as.matrix(Yct_b)
          if(!is.numeric(D_b)) mode(D_b) <- "numeric"
          if(!is.numeric(Y_b)) { warning(paste("Coercing replicate", b, "Y to numeric.")); mode(Y_b) <- "numeric" }
          if(!is.numeric(Yct_b)) { warning(paste("Coercing replicate", b, "Y.ct to numeric.")); mode(Yct_b) <- "numeric" }

          N_b <- ncol(D_b)
          expected_N_b <- N

          present_original_indices <- if (N > 0) 1:N else integer(0)
          replicate_col_to_original_idx <- if (N_b > 0) 1:N_b else integer(0) # Default mapping for bootstrap

          if(vartype == "jackknife") {
            expected_N_b <- N - 1
            if (N_b == N - 1 && n_replicates == N && N > 0) { # N > 0 for jackknife to make sense
              omitted_original_idx <- b
              present_original_indices <- (1:N)[-omitted_original_idx]
              # Map columns in this jackknife replicate D_b back to their original unit indices
              # Assuming D_b's columns are the N-1 present units, in their original relative order.
              replicate_col_to_original_idx <- present_original_indices
            } else if (N > 0) { # Only warn if N > 0, jackknife on N=0 is already handled.
              warning(paste("Replicate", b, ": Jackknife replicate has", N_b, "columns (expected N-1=", N-1, ") or n_replicates (", n_replicates, ") != N (", N, "). Skipping replicate."), call. = FALSE)
              next
            }
            # if N=0, expected_N_b = -1. N_b is likely 0. This path is complex if N=0.
            # However, N=0 is caught by an early return. So N > 0 here.
          }


          if (N > 0 && N_b != expected_N_b) { # Check for N > 0 before N_b != expected_N_b
            warning(paste("Replicate", b, ": For vartype '", vartype, "', expected ", expected_N_b, " columns in D_b, but got ", N_b, ". Skipping replicate."), call. = FALSE)
            next
          }


          if (nrow(D_b) != T_val || nrow(Y_b) != T_val || nrow(Yct_b) != T_val || (N_b > 0 && (ncol(Y_b) != N_b || ncol(Yct_b) != N_b))) {
            warning(paste("Replicate", b, ": Dimension mismatch with T_val. Skipping replicate."), call. = FALSE)
            next
          }

          if (N_b == 0 && N > 0) { next } # If N > 0, but this replicate has 0 units, skip. If N=0, N_b=0 is expected.

          Y_b[is.na(Yct_b)] <- NA


          # Get the column indices *in the current replicate D_b* that correspond to these target units
          replicate_cols_for_target_units <- which(colSums(D_b == 1, na.rm = TRUE) > 0)
          if (length(replicate_cols_for_target_units) == 0 && Ntr > 0) {
            next
          }


          for (t in seq.int(T_val)) {
            cols_to_avg <- integer(0)
            if (Ntr == 0) { # No treated units in D0, so effectively averaging "control" outcomes for "placebo" treated group
            } else if (!is.na(common_G_D0) && t < common_G_D0) { # Pre-treatment period for non-staggered
              cols_to_avg <- replicate_cols_for_target_units
            } else if (!is.na(common_G_D0) && t >= common_G_D0) { # Post-treatment period for non-staggered
              # Only average over units that are actually treated (D_b[t,col]==1) in this replicate at this time
              cols_to_avg <- D_b[t,] == 1
            } else if (is.na(common_G_D0) && Ntr > 0) { # Should not happen if is_staggered is FALSE and Ntr > 0
              warning(paste("Replicate", b, ": In non-staggered case with Ntr > 0, but common_G_D0 is NA. Skipping time", t, "for this replicate."), call. = FALSE)
              next # Skip this time period for this replicate
            }
            # If Ntr > 0 but replicate_cols_for_target_units is empty, cols_to_avg will be empty.

            if (length(cols_to_avg) > 0) {
              calendar_repl_tr[t, b] <- mean(Y_b[t, cols_to_avg], na.rm = TRUE)
              calendar_repl_cf[t, b] <- mean(Yct_b[t, cols_to_avg], na.rm = TRUE)
            }
            # If cols_to_avg is empty, the NA_real_ initialized in calendar_repl_tr/cf remains.
          }
        }
      }

      # MODIFICATION START: Keep all calendar times; rows with all NA replicates will result in NA averages.
      periods <- calendar_periods_labels # Use the potentially custom labels
      repl_tr <- calendar_repl_tr
      repl_cf <- calendar_repl_cf

      # rowMeans with na.rm=TRUE produces NaN if all replicates for a period are NA.
      # These NaNs will propagate to NAs in the final data.frame.
      theta_tr <- rowMeans(repl_tr, na.rm = TRUE)
      theta_cf <- rowMeans(repl_cf, na.rm = TRUE)
      # MODIFICATION END

    } # End Staggered/Non-Staggered block

    # --- Confidence Interval Calculation (Non-Parametric) ---
    num_periods_out <- length(periods)

    ci_tr95 <- matrix(NA_real_, nrow = num_periods_out, ncol = 2, dimnames = list(NULL, c("lower","upper")))
    ci_cf95 <- matrix(NA_real_, nrow = num_periods_out, ncol = 2, dimnames = list(NULL, c("lower","upper")))
    ci_tr90 <- matrix(NA_real_, nrow = num_periods_out, ncol = 2, dimnames = list(NULL, c("lower","upper")))
    ci_cf90 <- matrix(NA_real_, nrow = num_periods_out, ncol = 2, dimnames = list(NULL, c("lower","upper")))

    if (vartype == "bootstrap") {
      if (n_replicates > 0 && num_periods_out > 0) {
        if (nrow(repl_tr) != num_periods_out || nrow(repl_cf) != num_periods_out) {
          # This warning should ideally not be triggered if logic is correct upstream.
          warning("Bootstrap replicate matrix dimensions mismatch with number of output periods. Skipping CI.")
        } else {
          ci_tr95 <- basic_ci_alpha(theta_tr, repl_tr, alpha = 0.05)
          ci_cf95 <- basic_ci_alpha(theta_cf, repl_cf, alpha = 0.05)
          ci_tr90 <- basic_ci_alpha(theta_tr, repl_tr, alpha = 0.10)
          ci_cf90 <- basic_ci_alpha(theta_cf, repl_cf, alpha = 0.10)
        }
      } else { if (num_periods_out > 0) warning("Cannot calculate bootstrap CIs: No replicates remain or available.") } # num_periods_out > 0 condition added to warning

    } else if (vartype == "jackknife") {
      if (n_replicates > 1 && num_periods_out > 0) { # Jackknife needs n_replicates > 1 (i.e., N > 1)
        if (nrow(repl_tr) != num_periods_out || nrow(repl_cf) != num_periods_out) {
          warning("Jackknife replicate matrix dimensions mismatch with number of output periods. Skipping CI.")
        } else {
          jackknife_ci_alpha <- function(theta, replicates, alpha) {
            n_p <- length(theta) # number of periods
            n_repl <- ncol(replicates) # number of jackknife replicates (should be N)

            # Handle cases where theta might be NA (e.g. all replicates for that period were NA)
            # or where all replicates for a period are NA
            # dev_sq will propagate NAs if theta is NA or replicates are NA
            dev_sq <- (replicates - theta)^2 # (num_periods_out x n_repl) matrix

            # sum_dev_sq will be NA if theta was NA.
            # If theta is a number, but all replicates[p,] are NA, sum_dev_sq[p] will be 0 (due to na.rm=TRUE).
            sum_dev_sq   <- rowSums(dev_sq, na.rm = TRUE)
            valid_counts <- rowSums(!is.na(replicates)) # Count non-NA replicates per period

            jk_variance <- ((n_repl - 1) / n_repl) * sum_dev_sq
            # For periods where all replicates were NA (so theta is NA), jk_variance is NA.
            # For periods where theta is a number but all replicates were NA (sum_dev_sq=0), jk_variance is 0.
            # For periods where valid_counts <=1, variance is not well-defined or unstable.
            jk_variance[valid_counts <= 1] <- NA_real_

            jk_se <- sqrt(jk_variance)
            jk_se[!is.finite(jk_se)] <- NA_real_ # Catches NaN from sqrt(negative) or Inf

            z_crit <- qnorm(1 - alpha/2)

            # ok condition: theta must be non-NA, se must be non-NA, and se must be positive.
            # If se is 0 (e.g. from sum_dev_sq=0), CI would be [theta, theta], which is not informative.
            # Setting to NA is more appropriate for such cases.
            ok    <- !is.na(theta) & !is.na(jk_se) & jk_se > 0

            lower <- rep(NA_real_, n_p)
            upper <- rep(NA_real_, n_p)
            lower[ok] <- theta[ok] - z_crit * jk_se[ok]
            upper[ok] <- theta[ok] + z_crit * jk_se[ok]

            ci <- cbind(lower, upper)
            colnames(ci) <- c("lower", "upper")
            ci
          }
          ci_tr95 <- jackknife_ci_alpha(theta_tr, repl_tr, alpha = 0.05)
          ci_cf95 <- jackknife_ci_alpha(theta_cf, repl_cf, alpha = 0.05)
          ci_tr90 <- jackknife_ci_alpha(theta_tr, repl_tr, alpha = 0.10)
          ci_cf90 <- jackknife_ci_alpha(theta_cf, repl_cf, alpha = 0.10)
        }
      } else { if (num_periods_out > 0) warning("Cannot calculate jackknife CIs: Need >1 replicate (N > 1) or no periods remain.") } # num_periods_out > 0 condition added
    } else {
      stop("Internal error: Invalid 'vartype' for non-parametric CI calculation.")
    }

    # --- Assemble into out$Y.avg ---
    if (!exists("out")) out <- list()

    # Check consistency of component lengths before creating data.frame
    # This was already robust, but good to be aware with the changes.
    # num_periods_out is now full length of event_times_range or calendar_periods_labels.
    # theta_tr/cf will have this length, possibly with NaNs.
    # CIs are initialized to this size with NAs.

    # Ensure CI matrices are correctly dimensioned and numeric, even if all NA
    ensure_numeric_matrix_cols <- function(mat, expected_rows) {
      if (is.null(mat) || !is.matrix(mat) || ncol(mat) != 2 || nrow(mat) != expected_rows) {
        # If mat is not as expected (e.g. NULL from a failed CI calculation, or wrong dims), return a compliant NA matrix
        return(matrix(NA_real_, nrow = expected_rows, ncol = 2, dimnames = list(NULL, c("lower","upper"))))
      }
      if (!is.numeric(mat)) { # Coerce if not numeric (e.g. character matrix)
        mat_num <- matrix(as.numeric(mat), nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))
        mat <- mat_num
      }
      mat[is.nan(mat)] <- NA_real_ # Standardize NaN to NA
      mat
    }

    ci_tr95 <- ensure_numeric_matrix_cols(ci_tr95, num_periods_out)
    ci_cf95 <- ensure_numeric_matrix_cols(ci_cf95, num_periods_out)
    ci_tr90 <- ensure_numeric_matrix_cols(ci_tr90, num_periods_out)
    ci_cf90 <- ensure_numeric_matrix_cols(ci_cf90, num_periods_out)

    if (num_periods_out > 0) {
      # Convert theta_tr/cf (which might contain NaN from rowMeans) to numeric vectors.
      # data.frame conversion will turn NaNs into NAs.
      theta_tr_vec <- if(length(theta_tr) == num_periods_out) as.numeric(theta_tr) else rep(NA_real_, num_periods_out)
      theta_cf_vec <- if(length(theta_cf) == num_periods_out) as.numeric(theta_cf) else rep(NA_real_, num_periods_out)

      out$Y.avg <- data.frame(
        period         = periods, # This now contains all original periods
        treated        = theta_tr_vec, # Will be NA if all replicates were NA for this period
        counterfactual = theta_cf_vec, # Will be NA if all replicates were NA for this period
        lower.tr       = ci_tr95[,"lower"], upper.tr       = ci_tr95[,"upper"],
        lower90.tr     = ci_tr90[,"lower"], upper90.tr     = ci_tr90[,"upper"],
        lower.ct       = ci_cf95[,"lower"], upper.ct       = ci_cf95[,"upper"],
        lower90.ct     = ci_cf90[,"lower"], upper90.ct     = ci_cf90[,"upper"]
      )
    } else { # This case handles T_val=0 or N=0 (via early return) or if periods somehow becomes length 0
      out$Y.avg <- data.frame(period=numeric(0), treated=numeric(0), counterfactual=numeric(0),
                              lower.tr=numeric(0), upper.tr=numeric(0), lower90.tr=numeric(0), upper90.tr=numeric(0),
                              lower.ct=numeric(0), upper.ct=numeric(0), lower90.ct=numeric(0), upper90.ct=numeric(0))
    }
    attr(out$Y.avg, "period_type") <- period_type_attr

  } # when vartype == "parametric" the processing is done in plot.R
  ####################################
  ## Variance and CIs
  ####################################

  ## function to get two-sided p-values
  get.pvalue <- function(vec) {
    if (NaN %in% vec | NA %in% vec) {
      nan.pos <- is.nan(vec)
      na.pos <- is.na(vec)
      pos <- c(which(nan.pos), which(na.pos))
      vec.a <- vec[-pos]
      a <- sum(vec.a >= 0) / (length(vec) - sum(nan.pos | na.pos)) * 2
      b <- sum(vec.a <= 0) / (length(vec) - sum(nan.pos | na.pos)) * 2
    } else {
      a <- sum(vec >= 0) / length(vec) * 2
      b <- sum(vec <= 0) / length(vec) * 2
    }
    return(min(as.numeric(min(a, b)), 1))
  }

  ## ATT estimates
  if (vartype == "jackknife") {
    att.j <- jackknifed(att, att.boot, alpha, quantile.CI = quantile.CI)
    att.bound <- cbind(att + qnorm(alpha) * att.j$se, att + qnorm(1 - alpha) * att.j$se)
    colnames(att.bound) <- c("CI.lower", "CI.upper")
    rownames(att.bound) <- out$time

    est.att <- cbind(att, att.j$se, att.j$CI.l, att.j$CI.u, att.j$P, out$count)
    est.att90 <- cbind(att, att.j$se, att.bound, att.j$P, out$count)
    colnames(est.att) <- c(
      "ATT", "S.E.", "CI.lower", "CI.upper",
      "p.value", "count"
    )
    colnames(est.att90) <- c(
      "ATT", "S.E.", "CI.lower", "CI.upper",
      "p.value", "count"
    )
    rownames(est.att) <- out$time
    rownames(est.att90) <- out$time
    vcov.att <- att.j$vcov


    eff.calendar.j <- jackknifed(calendar.eff, calendar.eff.boot, alpha, quantile.CI = quantile.CI)
    est.eff.calendar <- cbind(calendar.eff, eff.calendar.j$se, eff.calendar.j$CI.l, eff.calendar.j$CI.u, eff.calendar.j$P, calendar.N)
    colnames(est.eff.calendar) <- c("ATT-calendar", "S.E.", "CI.lower", "CI.upper", "p.value", "count")

    eff.calendar.fit.j <- jackknifed(calendar.eff.fit, calendar.eff.fit.boot, alpha, quantile.CI = quantile.CI)
    est.eff.calendar.fit <- cbind(calendar.eff.fit, eff.calendar.fit.j$se, eff.calendar.fit.j$CI.l, eff.calendar.fit.j$CI.u, eff.calendar.fit.j$P, calendar.N)
    colnames(est.eff.calendar.fit) <- c("ATT-calendar Fitted", "S.E.", "CI.lower", "CI.upper", "p.value", "count")

    if (hasRevs == 1) {
      att.off.j <- jackknifed(att.off, att.off.boot, alpha, quantile.CI = quantile.CI)
      est.att.off <- cbind(att.off, att.off.j$se, att.off.j$CI.l, att.off.j$CI.u, att.off.j$P, out$count.off)
      colnames(est.att.off) <- c(
        "ATT.OFF", "S.E.", "CI.lower", "CI.upper",
        "p.value", "count"
      )
      rownames(est.att.off) <- out$time.off
      vcov.att.off <- att.off.j$vcov

      att.off.bound <- cbind(att.off + qnorm(alpha) * att.off.j$se, att.off + qnorm(1 - alpha) * att.off.j$se)
      colnames(att.off.bound) <- c("CI.lower", "CI.upper")
      rownames(att.off.bound) <- out$time.off
    }

    if (!is.null(T.on.carry)) {
      carry.att.j <- jackknifed(carry.att, carry.att.boot, alpha, quantile.CI = quantile.CI)
      est.carry.att <- cbind(
        carry.att, carry.att.j$se,
        carry.att.j$CI.l, carry.att.j$CI.u, carry.att.j$P
      )

      colnames(est.carry.att) <- c(
        "ATT", "S.E.", "CI.lower", "CI.upper",
        "p.value"
      )
      rownames(est.carry.att) <- carry.time
    }

    if (!is.null(balance.period)) {
      balance.att.j <- jackknifed(balance.att, balance.att.boot, alpha, quantile.CI = quantile.CI)
      est.balance.att <- cbind(balance.att, balance.att.j$se, balance.att.j$CI.l, balance.att.j$CI.u, balance.att.j$P, out$balance.count)
      colnames(est.balance.att) <- c(
        "ATT", "S.E.", "CI.lower", "CI.upper",
        "p.value", "count"
      )
      rownames(est.balance.att) <- out$balance.time
      vcov.balance.att <- balance.att.j$vcov


      balance.avg.att.j <- jackknifed(balance.avg.att, balance.avg.att.boot, alpha, quantile.CI = quantile.CI)
      est.balance.avg <- t(as.matrix(c(balance.avg.att, balance.avg.att.j$se, balance.avg.att.j$CI.l, balance.avg.att.j$CI.u, balance.avg.att.j$P)))
      colnames(est.balance.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")

      balance.att.bound <- cbind(balance.att + qnorm(alpha) * balance.att.j$se, balance.att + qnorm(1 - alpha) * balance.att.j$se)
      colnames(balance.att.bound) <- c("CI.lower", "CI.upper")
      rownames(balance.att.bound) <- out$balance.time

      if (!is.null(placebo.period) & placeboTest == TRUE) {
        balance.att.placebo.j <- jackknifed(balance.att.placebo, balance.att.placebo.boot, alpha, quantile.CI = quantile.CI)
        est.balance.placebo <- t(as.matrix(c(balance.att.placebo, balance.att.placebo.j$se, balance.att.placebo.j$CI.l, balance.att.placebo.j$CI.u, balance.att.placebo.j$P)))
        colnames(est.balance.placebo) <- c("ATT.placebo", "S.E.", "CI.lower", "CI.upper", "p.value")
      }
    }

    if (!is.null(W)) {
      att.avg.W.j <- jackknifed(att.avg.W, att.avg.W.boot, alpha, quantile.CI = quantile.CI)
      est.avg.W <- t(as.matrix(c(att.avg.W, att.avg.W.j$se, att.avg.W.j$CI.l, att.avg.W.j$CI.u, att.avg.W.j$P)))
      colnames(est.avg.W) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")

      att.on.W.j <- jackknifed(att.on.W, att.on.W.boot, alpha, quantile.CI = quantile.CI)
      est.att.W <- cbind(att.on.W, att.on.W.j$se, att.on.W.j$CI.l, att.on.W.j$CI.u, att.on.W.j$P, count.on.W)
      colnames(est.att.W) <- c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value", "count")
      rownames(est.att.W) <- time.on.W

      vcov.att.W <- att.on.W.j$vcov

      att.W.bound <- cbind(att.on.W + qnorm(alpha) * att.on.W.j$se, att.on.W + qnorm(1 - alpha) * att.on.W.j$se)
      colnames(att.W.bound) <- c("CI.lower", "CI.upper")
      rownames(att.W.bound) <- time.on.W

      if (!is.null(placebo.period) & placeboTest == TRUE) {
        att.placebo.W.j <- jackknifed(att.placebo.W, att.placebo.W.boot, alpha, quantile.CI = quantile.CI)
        est.placebo.W <- t(as.matrix(c(att.placebo.W, att.placebo.W.j$se, att.placebo.W.j$CI.l, att.placebo.W.j$CI.u, att.placebo.W.j$P)))
        colnames(est.placebo.W) <- c("ATT.placebo", "S.E.", "CI.lower", "CI.upper", "p.value")
      }
      if (hasRevs == 1) {
        att.off.W.j <- jackknifed(att.off.W, att.off.W.boot, alpha, quantile.CI = quantile.CI)
        est.att.off.W <- cbind(att.off.W, att.off.W.j$se, att.off.W.j$CI.l, att.off.W.j$CI.u, att.off.W.j$P, count.off.W)
        colnames(est.att.off.W) <- c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value", "count")
        rownames(est.att.off.W) <- time.off.W
        vcov.att.off.W <- att.off.W.j$vcov
        att.off.W.bound <- cbind(att.off.W + qnorm(alpha) * att.off.W.j$se, att.off.W + qnorm(1 - alpha) * att.off.W.j$se)
        colnames(att.off.W.bound) <- c("CI.lower", "CI.upper")
        rownames(att.off.W.bound) <- out$time.off

        if (!is.null(carryover.period) & carryoverTest == TRUE) {
          att.carryover.W.j <- jackknifed(att.carryover.W, att.carryover.W.boot, alpha, quantile.CI = quantile.CI)
          est.carryover.W <- t(as.matrix(c(att.carryover.W, att.carryover.W.j$se, att.carryover.W.j$CI.l, att.carryover.W.j$CI.u, att.carryover.W.j$P)))
          colnames(est.carryover.W) <- c("ATT.carryover", "S.E.", "CI.lower", "CI.upper", "p.value")
        }
      }
    }

    ## average (over time) ATT
    att.avg.j <- jackknifed(att.avg, att.avg.boot, alpha, quantile.CI = quantile.CI)
    est.avg <- t(as.matrix(c(att.avg, att.avg.j$se, att.avg.j$CI.l, att.avg.j$CI.u, att.avg.j$P)))
    colnames(est.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")

    att.avg.unit.j <- jackknifed(att.avg.unit, att.avg.unit.boot, alpha, quantile.CI = quantile.CI)
    est.avg.unit <- t(as.matrix(c(att.avg.unit, att.avg.unit.j$se, att.avg.unit.j$CI.l, att.avg.unit.j$CI.u, att.avg.unit.j$P)))
    colnames(est.avg.unit) <- c("ATT.avg.unit", "S.E.", "CI.lower", "CI.upper", "p.value")

    ## regression coefficents
    if (p > 0) {
      beta.j <- jackknifed(beta, beta.boot, alpha, quantile.CI = quantile.CI)
      est.beta <- cbind(beta, beta.j$se, beta.j$CI.l, beta.j$CI.u, beta.j$P)
      colnames(est.beta) <- c("beta", "S.E.", "CI.lower", "CI.upper", "p.value")

      if (binary == TRUE) {
        marginal.j <- jackknifed(out$marginal, marginal.boot, alpha, quantile.CI = quantile.CI)
        est.marginal <- cbind(out$marginal, marginal.j$se, marginal.j$CI.l, marginal.j$CI.u, marginal.j$P)
        colnames(est.marginal) <- c("marginal", "S.E.", "CI.lower", "CI.upper", "p.value")
      }
    }

    ## placebo test
    if (!is.null(placebo.period) & placeboTest == TRUE) {
      att.placebo <- out$att.placebo
      att.placebo.j <- jackknifed(att.placebo, att.placebo.boot, alpha, quantile.CI = quantile.CI)
      att.placebo.bound <- c(
        att.placebo + qnorm(alpha) * att.placebo.j$se,
        att.placebo + qnorm(1 - alpha) * att.placebo.j$se
      )
      est.placebo <- t(as.matrix(c(
        att.placebo, att.placebo.j$se,
        att.placebo.j$CI.l, att.placebo.j$CI.u,
        att.placebo.j$P,
        att.placebo.bound
      )))
      colnames(est.placebo) <- c(
        "ATT.placebo", "S.E.",
        "CI.lower", "CI.upper",
        "p.value", "CI.lower(90%)", "CI.upper(90%)"
      )
    }

    ## carryover test
    if (!is.null(carryover.period) & carryoverTest == TRUE) {
      att.carryover <- out$att.carryover
      att.carryover.j <- jackknifed(att.carryover, att.carryover.boot, alpha, quantile.CI = quantile.CI)
      att.carryover.bound <- c(
        att.carryover + qnorm(alpha) * att.carryover.j$se,
        att.carryover + qnorm(1 - alpha) * att.carryover.j$se
      )

      est.carryover <- t(as.matrix(c(
        att.carryover, att.carryover.j$se,
        att.carryover.j$CI.l, att.carryover.j$CI.u,
        att.carryover.j$P,
        att.carryover.bound
      )))
      colnames(est.carryover) <- c(
        "ATT.carryover", "S.E.",
        "CI.lower", "CI.upper",
        "p.value", "CI.lower(90%)", "CI.upper(90%)"
      )
    }

    ## cohort effect
    est.group.out <- NULL
    if (!is.null(group)) {
      group.att.j <- jackknifed(group.att, group.att.boot, alpha, quantile.CI = quantile.CI)
      est.group.att <- cbind(group.att, group.att.j$se, group.att.j$CI.l, group.att.j$CI.u, group.att.j$P)
      colnames(est.group.att) <- c(
        "ATT", "S.E.", "CI.lower", "CI.upper",
        "p.value"
      )

      est.group.out <- list()
      for (sub.name in group.output.name) {
        subgroup.atts <- group.output.origin[[sub.name]]$att.on
        subgroup.atts.boot <- group.atts.boot[[sub.name]]
        subgroup.est.att <- NULL
        subgroup.att.bound <- NULL

        if (dim(subgroup.atts.boot)[1] > 0) {
          subgroup.att.j <- jackknifed(subgroup.atts, subgroup.atts.boot, alpha, quantile.CI = quantile.CI)
          subgroup.est.att <- cbind(
            subgroup.atts, subgroup.att.j$se, subgroup.att.j$CI.l,
            subgroup.att.j$CI.u, subgroup.att.j$P,
            group.output.origin[[sub.name]]$count.on
          )
          colnames(subgroup.est.att) <- c(
            "ATT", "S.E.", "CI.lower", "CI.upper",
            "p.value", "count"
          )
          rownames(subgroup.est.att) <- group.output.origin[[sub.name]]$time.on

          subgroup.att.bound <- cbind(
            subgroup.atts + qnorm(alpha) * subgroup.att.j$se,
            subgroup.atts + qnorm(1 - alpha) * subgroup.att.j$se
          )
          colnames(subgroup.att.bound) <- c("CI.lower", "CI.upper")
          rownames(subgroup.att.bound) <- group.output.origin[[sub.name]]$time.on
        }

        subgroup.est.att.off <- NULL
        subgroup.att.off.bound <- NULL
        if (hasRevs == 1) {
          subgroup.atts.off <- group.output.origin[[sub.name]]$att.off
          subgroup.atts.off.boot <- group.atts.off.boot[[sub.name]]
          if (dim(subgroup.atts.off.boot)[1] > 0) {
            subgroup.att.off.j <- jackknifed(subgroup.atts.off, subgroup.atts.off.boot, alpha, quantile.CI = quantile.CI)
            subgroup.est.att.off <- cbind(
              subgroup.atts.off, subgroup.att.off.j$se, subgroup.att.off.j$CI.l,
              subgroup.att.off.j$CI.u, subgroup.att.off.j$P,
              group.output.origin[[sub.name]]$count.off
            )
            colnames(subgroup.est.att.off) <- c(
              "ATT", "S.E.", "CI.lower", "CI.upper",
              "p.value", "count"
            )
            rownames(subgroup.est.att.off) <- group.output.origin[[sub.name]]$time.off

            subgroup.att.off.bound <- cbind(
              subgroup.atts.off + qnorm(alpha) * subgroup.att.off.j$se,
              subgroup.atts.off + qnorm(1 - alpha) * subgroup.att.off.j$se
            )
            colnames(subgroup.att.off.bound) <- c("CI.lower", "CI.upper")
            rownames(subgroup.att.off.bound) <- group.output.origin[[sub.name]]$time.off
          }
        }

        subgroup.est.placebo <- NULL
        if (placeboTest) {
          subgroup.att.placebo <- group.output.origin[[sub.name]]$att.placebo
          if (length(subgroup.att.placebo) > 0) {
            subgroup.att.placebo.j <- jackknifed(subgroup.att.placebo, group.att.placebo.boot[[sub.name]], alpha, quantile.CI = quantile.CI)
            att.placebo.bound <- c(
              subgroup.att.placebo + qnorm(alpha) * subgroup.att.placebo.j$se,
              subgroup.att.placebo + qnorm(1 - alpha) * subgroup.att.placebo.j$se
            )

            subgroup.est.placebo <- t(as.matrix(c(
              subgroup.att.placebo,
              subgroup.att.placebo.j$se,
              subgroup.att.placebo.j$CI.l,
              subgroup.att.placebo.j$CI.u,
              subgroup.att.placebo.j$P,
              att.placebo.bound
            )))
            colnames(subgroup.est.placebo) <- c(
              "ATT.placebo", "S.E.",
              "CI.lower", "CI.upper", "p.value",
              "CI.lower(90%)", "CI.upper(90%)"
            )
          }
        }

        subgroup.est.carryover <- NULL
        if (carryoverTest) {
          subgroup.att.carryover <- group.output.origin[[sub.name]]$att.carryover
          if (length(subgroup.att.carryover) > 0) {
            subgroup.att.carryover.j <- jackknifed(subgroup.att.carryover, group.att.carryover.boot[[sub.name]], alpha, quantile.CI = quantile.CI)
            att.carryover.bound <- c(
              subgroup.att.carryover + qnorm(alpha) * subgroup.att.carryover.j$se,
              subgroup.att.carryover + qnorm(1 - alpha) * subgroup.att.carryover.j$se
            )

            subgroup.est.carryover <- t(as.matrix(c(
              subgroup.att.carryover,
              subgroup.att.carryover.j$se,
              subgroup.att.carryover.j$CI.l,
              subgroup.att.carryover.j$CI.u,
              subgroup.att.carryover.j$P,
              att.carryover.bound
            )))
            colnames(subgroup.est.carryover) <- c(
              "ATT.carryover", "S.E.",
              "CI.lower", "CI.upper", "p.value",
              "CI.lower(90%)", "CI.upper(90%)"
            )
          }
        }

        est.group.out[[sub.name]] <- list(
          att.on = subgroup.est.att,
          att.on.bound = subgroup.att.bound,
          att.on.boot = group.atts.boot[[sub.name]],
          att.off = subgroup.est.att.off,
          att.off.bound = subgroup.att.off.bound,
          att.off.boot = group.atts.off.boot[[sub.name]],
          att.placebo = subgroup.est.placebo,
          att.carryover = subgroup.est.carryover
        )
      }
    }
  } else {
    se.att <- apply(att.boot, 1, function(vec) sd(vec, na.rm = TRUE))
    if (quantile.CI == FALSE) {
      CI.att <- cbind(att - se.att * qnorm(1 - alpha / 2), att + se.att * qnorm(1 - alpha / 2)) # normal approximation
      pvalue.att <- (1 - pnorm(abs(att / se.att))) * 2
    } else {
      CI.att <- t(apply(att.boot, 1, function(vec) {
        2 * att[which.max(!is.na(vec))] - quantile(vec, c(1 - alpha / 2, alpha / 2), na.rm = TRUE)
      }))
      pvalue.att <- apply(att.boot, 1, get.pvalue)
    }

    #vcov.att <- cov(t(att.boot), use = "pairwise.complete.obs")
    vcov.att <- tryCatch(
      {
        cov(t(att.boot), use = "pairwise.complete.obs")
      },
      error = function(e) {
        NA
      }
    )
    # for equivalence test
    if (quantile.CI == FALSE) {
      att.bound <- cbind(att - se.att * qnorm(1 - alpha), att + se.att * qnorm(1 - alpha)) # one-sided
    } else {
      att.bound <- t(apply(att.boot, 1, function(vec) quantile(vec, c(alpha, 1 - alpha), na.rm = TRUE)))
    }

    colnames(att.bound) <- c("CI.lower", "CI.upper")
    rownames(att.bound) <- out$time
    est.att <- cbind(att, se.att, CI.att, pvalue.att, out$count)
    est.att90 <- cbind(att, se.att, att.bound, pvalue.att, out$count)
    colnames(est.att) <- c(
      "ATT", "S.E.", "CI.lower", "CI.upper",
      "p.value", "count"
    )
    colnames(est.att90) <- c(
      "ATT", "S.E.", "CI.lower", "CI.upper",
      "p.value", "count"
    )
    rownames(est.att) <- out$time
    rownames(est.att90) <- out$time




    if (hasRevs == 1) {
      se.att.off <- apply(att.off.boot, 1, function(vec) sd(vec, na.rm = TRUE))
      if (quantile.CI == FALSE) {
        CI.att.off <- cbind(att.off - se.att.off * qnorm(1 - alpha / 2), att.off + se.att.off * qnorm(1 - alpha / 2))
        pvalue.att.off <- (1 - pnorm(abs(att.off / se.att.off))) * 2
      } else {
        CI.att.off <- t(apply(att.off.boot, 1, function(vec) quantile(vec, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)))
        pvalue.att.off <- apply(att.off.boot, 1, get.pvalue)
      }

      #vcov.att.off <- cov(t(att.off.boot), use = "pairwise.complete.obs")
      vcov.att.off <- tryCatch(
        {
          cov(t(att.off.boot), use = "pairwise.complete.obs")
        },
        error = function(e) {
          NA
        }
      )


      est.att.off <- cbind(att.off, se.att.off, CI.att.off, pvalue.att.off, out$count.off)
      colnames(est.att.off) <- c(
        "ATT.OFF", "S.E.", "CI.lower", "CI.upper",
        "p.value", "count.off"
      )
      rownames(est.att.off) <- out$time.off
      # T0.off.l <- sum(out$time.off > 0)
      # norm.att.off.sq <- (att.off/se.att.off)^2
      # T0.off.p <- 1 - pchisq(sum(norm.att.off.sq[(length(out$time.off) - T0.off.l + 1):length(out$time.off)]), df = T0.off.l)
      if (quantile.CI == FALSE) {
        att.off.bound <- cbind(att.off - se.att.off * qnorm(1 - alpha), att.off + se.att.off * qnorm(1 - alpha))
      } else {
        att.off.bound <- t(apply(att.off.boot, 1, function(vec) quantile(vec, c(alpha, 1 - alpha), na.rm = TRUE)))
      }

      colnames(att.off.bound) <- c("CI.lower", "CI.upper")
      rownames(att.off.bound) <- out$time.off
    }

    if (!is.null(T.on.carry)) {
      se.carry.att <- apply(carry.att.boot, 1, function(vec) sd(vec, na.rm = TRUE))
      if (quantile.CI == FALSE) {
        CI.carry.att <- cbind(
          carry.att - se.carry.att * qnorm(1 - alpha / 2),
          carry.att + se.carry.att * qnorm(1 - alpha / 2)
        ) # normal approximation
        pvalue.carry.att <- (1 - pnorm(abs(carry.att / se.carry.att))) * 2
      } else {
        CI.carry.att <- t(apply(carry.att.boot, 1, function(vec) quantile(vec, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)))
        pvalue.carry.att <- apply(carry.att.boot, 1, get.pvalue)
      }


      est.carry.att <- cbind(
        carry.att, se.carry.att,
        CI.carry.att, pvalue.carry.att
      )

      colnames(est.carry.att) <- c(
        "ATT", "S.E.", "CI.lower", "CI.upper",
        "p.value"
      )
      rownames(est.carry.att) <- carry.time
    }

    if (!is.null(balance.period)) {
      se.balance.att <- apply(balance.att.boot, 1, function(vec) sd(vec, na.rm = TRUE))
      if (quantile.CI == FALSE) {
        CI.balance.att <- cbind(
          balance.att - se.balance.att * qnorm(1 - alpha / 2),
          balance.att + se.balance.att * qnorm(1 - alpha / 2)
        )
        pvalue.balance.att <- (1 - pnorm(abs(balance.att / se.balance.att))) * 2
      } else {
        CI.balance.att <- t(apply(balance.att.boot, 1, function(vec) quantile(vec, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)))
        pvalue.balance.att <- apply(balance.att.boot, 1, get.pvalue)
      }

      #vcov.balance.att <- cov(t(balance.att.boot), use = "pairwise.complete.obs")

      vcov.balance.att <- tryCatch(
        {
          cov(t(balance.att.boot), use = "pairwise.complete.obs")
        },
        error = function(e) {
          NA
        }
      )

      est.balance.att <- cbind(
        balance.att, se.balance.att, CI.balance.att,
        pvalue.balance.att, out$balance.count
      )
      colnames(est.balance.att) <- c(
        "ATT", "S.E.", "CI.lower", "CI.upper",
        "p.value", "count"
      )
      rownames(est.balance.att) <- out$balance.time

      se.balance.avg.att <- sd(balance.avg.att.boot, na.rm = TRUE)
      if (quantile.CI == FALSE) {
        CI.balance.avg.att <- c(
          balance.avg.att - se.balance.avg.att * qnorm(1 - alpha / 2),
          balance.avg.att + se.balance.avg.att * qnorm(1 - alpha / 2)
        )
        p.balance.avg.att <- (1 - pnorm(abs(balance.avg.att / se.balance.avg.att))) * 2
      } else {
        CI.balance.avg.att <- quantile(balance.avg.att.boot, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
        p.balance.avg.att <- get.pvalue(balance.avg.att.boot)
      }

      est.balance.avg <- t(as.matrix(c(balance.avg.att, se.balance.avg.att, CI.balance.avg.att, p.balance.avg.att)))
      colnames(est.balance.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")
      if (quantile.CI == FALSE) {
        balance.att.bound <- cbind(
          balance.att - se.balance.att * qnorm(1 - alpha),
          balance.att + se.balance.att * qnorm(1 - alpha)
        )
      } else {
        balance.att.bound <- t(apply(balance.att.boot, 1, function(vec) quantile(vec, c(alpha, 1 - alpha), na.rm = TRUE)))
      }

      colnames(balance.att.bound) <- c("CI.lower", "CI.upper")
      rownames(balance.att.bound) <- out$balance.time

      if (!is.null(placebo.period) & placeboTest == TRUE) {
        balance.att.placebo <- out$balance.att.placebo
        balance.se.placebo <- sd(balance.att.placebo.boot, na.rm = TRUE)
        if (quantile.CI == FALSE) {
          balance.CI.placebo <- c(
            balance.att.placebo - balance.se.placebo * qnorm(1 - alpha / 2),
            balance.att.placebo + balance.se.placebo * qnorm(1 - alpha / 2)
          )
          balance.CI.placebo.bound <- c(
            balance.att.placebo - balance.se.placebo * qnorm(1 - alpha),
            balance.att.placebo + balance.se.placebo * qnorm(1 - alpha)
          )
          balance.pvalue.placebo <- (1 - pnorm(abs(balance.att.placebo / balance.se.placebo))) * 2
        } else {
          balance.CI.placebo <- quantile(balance.att.placebo.boot, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
          balance.CI.placebo.bound <- quantile(balance.att.placebo.boot, c(alpha, 1 - alpha), na.rm = TRUE)
          balance.pvalue.placebo <- get.pvalue(balance.att.placebo.boot)
        }


        est.balance.placebo <- t(as.matrix(c(
          balance.att.placebo,
          balance.se.placebo,
          balance.CI.placebo,
          balance.pvalue.placebo,
          balance.CI.placebo.bound
        )))
        colnames(est.balance.placebo) <- c(
          "ATT.placebo", "S.E.",
          "CI.lower", "CI.upper", "p.value",
          "CI.lower(90%)", "CI.upper(90%)"
        )
      }
    }

    if (!is.null(W)) {
      # att.avg.W.boot
      se.att.avg.W <- sd(att.avg.W.boot, na.rm = TRUE)
      if (quantile.CI == FALSE) {
        CI.att.avg.W <- c(
          att.avg.W - se.att.avg.W * qnorm(1 - alpha / 2),
          att.avg.W + se.att.avg.W * qnorm(1 - alpha / 2)
        )
        p.att.avg.W <- (1 - pnorm(abs(att.avg.W / se.att.avg.W))) * 2
      } else {
        CI.att.avg.W <- quantile(att.avg.W.boot, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
        p.att.avg.W <- get.pvalue(att.avg.W.boot)
      }

      est.avg.W <- t(as.matrix(c(att.avg.W, se.att.avg.W, CI.att.avg.W, p.att.avg.W)))
      colnames(est.avg.W) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")

      # att.on.W.boot
      se.att.W <- apply(att.on.W.boot, 1, function(vec) sd(vec, na.rm = TRUE))
      if (quantile.CI == FALSE) {
        CI.att.W <- cbind(
          att.on.W - se.att.W * qnorm(1 - alpha / 2),
          att.on.W + se.att.W * qnorm(1 - alpha / 2)
        )
        att.W.bound <- cbind(
          att.on.W - se.att.W * qnorm(1 - alpha),
          att.on.W + se.att.W * qnorm(1 - alpha)
        )
        pvalue.att.W <- (1 - pnorm(abs(att.on.W / se.att.W))) * 2
      } else {
        CI.att.W <- t(apply(att.on.W.boot, 1, function(vec) quantile(vec, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)))
        att.W.bound <- t(apply(att.on.W.boot, 1, function(vec) quantile(vec, c(alpha, 1 - alpha), na.rm = TRUE)))
        pvalue.att.W <- apply(att.on.W.boot, 1, get.pvalue)
      }

      #vcov.att.W <- cov(t(att.on.W.boot), use = "pairwise.complete.obs")
      vcov.att.W <- tryCatch(
        {
          cov(t(att.on.W.boot), use = "pairwise.complete.obs")
        },
        error = function(e) {
          NA
        }
      )

      est.att.W <- cbind(
        att.on.W, se.att.W, CI.att.W,
        pvalue.att.W, count.on.W
      )
      colnames(est.att.W) <- c(
        "ATT", "S.E.", "CI.lower", "CI.upper",
        "p.value", "count"
      )
      rownames(est.att.W) <- time.on.W


      colnames(att.W.bound) <- c("CI.lower", "CI.upper")
      rownames(att.W.bound) <- time.on.W

      if (!is.null(placebo.period) & placeboTest == TRUE) {
        # att.placebo.W.boot
        se.placebo.W <- sd(att.placebo.W.boot, na.rm = TRUE)
        if (quantile.CI == FALSE) {
          CI.placebo.W <- c(
            att.placebo.W - se.placebo.W * qnorm(1 - alpha / 2),
            att.placebo.W + se.placebo.W * qnorm(1 - alpha / 2)
          )
          CI.placebo.bound.W <- c(
            att.placebo.W - se.placebo.W * qnorm(1 - alpha),
            att.placebo.W + se.placebo.W * qnorm(1 - alpha)
          )
          pvalue.placebo.w <- (1 - pnorm(abs(att.placebo.W / se.placebo.W))) * 2
        } else {
          CI.placebo.W <- quantile(att.placebo.W.boot, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)

          CI.placebo.bound.W <- quantile(att.placebo.W.boot, c(alpha, 1 - alpha), na.rm = TRUE)

          pvalue.placebo.w <- get.pvalue(att.placebo.W.boot)
        }

        est.placebo.W <- t(as.matrix(c(
          att.placebo.W,
          se.placebo.W,
          CI.placebo.W,
          pvalue.placebo.w,
          CI.placebo.bound.W
        )))
        colnames(est.placebo.W) <- c(
          "ATT.placebo", "S.E.",
          "CI.lower", "CI.upper", "p.value",
          "CI.lower(90%)", "CI.upper(90%)"
        )
      }
      if (hasRevs == 1) {
        # att.off.W.boot
        se.att.off.W <- apply(att.off.W.boot, 1, function(vec) sd(vec, na.rm = TRUE))
        if (quantile.CI == FALSE) {
          CI.att.off.W <- cbind(
            att.off.W - se.att.off.W * qnorm(1 - alpha / 2),
            att.off.W + se.att.off.W * qnorm(1 - alpha / 2)
          )
          att.off.W.bound <- cbind(
            att.off.W - se.att.off.W * qnorm(1 - alpha),
            att.off.W + se.att.off.W * qnorm(1 - alpha)
          )
          pvalue.att.off.W <- (1 - pnorm(abs(att.off.W / se.att.off.W))) * 2
        } else {
          CI.att.off.W <- t(apply(att.off.W.boot, 1, function(vec) quantile(vec, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)))
          att.off.W.bound <- t(apply(att.off.W.boot, 1, function(vec) quantile(vec, c(alpha, 1 - alpha), na.rm = TRUE)))
          pvalue.att.off.W <- apply(att.off.W.boot, 1, get.pvalue)
        }

        #vcov.att.off.W <- cov(t(att.off.W.boot), use = "pairwise.complete.obs")
        vcov.att.off.W <- tryCatch(
          {
            cov(t(att.off.W.boot), use = "pairwise.complete.obs")
          },
          error = function(e) {
            NA
          }
        )

        est.att.off.W <- cbind(
          att.off.W, se.att.off.W, CI.att.off.W,
          pvalue.att.off.W, count.off.W
        )
        colnames(est.att.off.W) <- c(
          "ATT", "S.E.", "CI.lower", "CI.upper",
          "p.value", "count"
        )
        rownames(est.att.off.W) <- time.off.W
        colnames(att.off.W.bound) <- c("CI.lower", "CI.upper")
        rownames(att.off.W.bound) <- time.off.W

        if (!is.null(carryover.period) & carryoverTest == TRUE) {
          # att.carryover.W.boot
          se.carryover.W <- sd(att.carryover.W.boot, na.rm = TRUE)
          if (quantile.CI == FALSE) {
            CI.carryover.W <- c(
              att.carryover.W - se.carryover.W * qnorm(1 - alpha / 2),
              att.carryover.W + se.carryover.W * qnorm(1 - alpha / 2)
            )
            CI.carryover.bound.W <- c(
              att.carryover.W - se.carryover.W * qnorm(1 - alpha),
              att.carryover.W + se.carryover.W * qnorm(1 - alpha)
            )
            pvalue.carryover.w <- (1 - pnorm(abs(att.carryover.W / se.carryover.W))) * 2
          } else {
            CI.carryover.W <- quantile(att.carryover.W.boot, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
            CI.carryover.bound.W <- quantile(att.carryover.W.boot, c(alpha, 1 - alpha), na.rm = TRUE)
            pvalue.carryover.w <- get.pvalue(att.carryover.W.boot)
          }

          est.carryover.W <- t(as.matrix(c(
            att.carryover.W,
            se.carryover.W,
            CI.carryover.W,
            pvalue.carryover.w,
            CI.carryover.bound.W
          )))
          colnames(est.carryover.W) <- c(
            "ATT.carryover", "S.E.",
            "CI.lower", "CI.upper", "p.value",
            "CI.lower(90%)", "CI.upper(90%)"
          )
        }
      }
    }

    ## average (over time) ATT
    se.avg <- sd(att.avg.boot, na.rm = TRUE)
    if (quantile.CI == FALSE) {
      CI.avg <- c(att.avg - se.avg * qnorm(1 - alpha / 2), att.avg + se.avg * qnorm(1 - alpha / 2))
      pvalue.avg <- (1 - pnorm(abs(att.avg / se.avg))) * 2
    } else {
      CI.avg <- quantile(att.avg.boot, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
      pvalue.avg <- get.pvalue(att.avg.boot)
    }

    est.avg <- t(as.matrix(c(att.avg, se.avg, CI.avg, pvalue.avg)))
    colnames(est.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")

    se.avg.unit <- sd(att.avg.unit.boot, na.rm = TRUE)
    if (quantile.CI == FALSE) {
      CI.avg.unit <- c(
        att.avg.unit - se.avg.unit * qnorm(1 - alpha / 2),
        att.avg.unit + se.avg.unit * qnorm(1 - alpha / 2)
      )
      pvalue.avg.unit <- (1 - pnorm(abs(att.avg.unit / se.avg.unit))) * 2
    } else {
      CI.avg.unit <- quantile(att.avg.unit.boot, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
      pvalue.avg.unit <- get.pvalue(att.avg.unit.boot)
    }

    est.avg.unit <- t(as.matrix(c(att.avg.unit, se.avg.unit, CI.avg.unit, pvalue.avg.unit)))
    colnames(est.avg.unit) <- c("ATT.avg.unit", "S.E.", "CI.lower", "CI.upper", "p.value")


    se.eff.calendar <- apply(calendar.eff.boot, 1, function(vec) sd(vec, na.rm = TRUE))
    if (quantile.CI == FALSE) {
      CI.eff.calendar <- cbind(calendar.eff - se.eff.calendar * qnorm(1 - alpha / 2), calendar.eff + se.eff.calendar * qnorm(1 - alpha / 2))
      pvalue.eff.calendar <- (1 - pnorm(abs(calendar.eff / se.eff.calendar))) * 2
    } else {
      CI.eff.calendar <- t(apply(calendar.eff.boot, 1, function(vec) quantile(vec, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)))
      pvalue.eff.calendar <- apply(calendar.eff.boot, 1, get.pvalue)
    }
    est.eff.calendar <- cbind(calendar.eff, se.eff.calendar, CI.eff.calendar, pvalue.eff.calendar, calendar.N)
    colnames(est.eff.calendar) <- c("ATT-calendar", "S.E.", "CI.lower", "CI.upper", "p.value", "count")

    se.eff.calendar.fit <- apply(calendar.eff.fit.boot, 1, function(vec) sd(vec, na.rm = TRUE))
    if (quantile.CI == FALSE) {
      CI.eff.calendar.fit <- cbind(calendar.eff.fit - se.eff.calendar.fit * qnorm(1 - alpha / 2), calendar.eff.fit + se.eff.calendar.fit * qnorm(1 - alpha / 2))
      pvalue.eff.calendar.fit <- (1 - pnorm(abs(calendar.eff.fit / se.eff.calendar.fit))) * 2
    } else {
      CI.eff.calendar.fit <- t(apply(calendar.eff.fit.boot, 1, function(vec) quantile(vec, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)))
      pvalue.eff.calendar.fit <- apply(calendar.eff.fit.boot, 1, get.pvalue)
    }
    est.eff.calendar.fit <- cbind(calendar.eff.fit, se.eff.calendar.fit, CI.eff.calendar.fit, pvalue.eff.calendar.fit, calendar.N)
    colnames(est.eff.calendar.fit) <- c("ATT-calendar Fitted", "S.E.", "CI.lower", "CI.upper", "p.value", "count")

    ## regression coefficents
    if (p > 0) {
      se.beta <- apply(beta.boot, 1, function(vec) sd(vec, na.rm = TRUE))
      if (quantile.CI == FALSE) {
        CI.beta <- cbind(c(beta) - se.beta * qnorm(1 - alpha / 2), c(beta) + se.beta * qnorm(1 - alpha / 2))
        pvalue.beta <- (1 - pnorm(abs(beta / se.beta))) * 2
      } else {
        CI.beta <- t(apply(beta.boot, 1, function(vec) quantile(vec, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)))
        pvalue.beta <- apply(beta.boot, 1, get.pvalue)
      }
      est.beta <- cbind(c(beta), se.beta, CI.beta, pvalue.beta)
      colnames(est.beta) <- c("Coef", "S.E.", "CI.lower", "CI.upper", "p.value")

      if (binary == TRUE) {
        out$marginal[na.pos] <- NA
        se.marginal <- apply(marginal.boot, 1, function(vec) sd(vec, na.rm = TRUE))
        CI.marginal <- cbind(
          c(out$marginal) - se.marginal * qnorm(1 - alpha / 2),
          c(out$marginal) + se.marginal * qnorm(1 - alpha / 2)
        )
        pvalue.marginal <- (1 - pnorm(abs(out$marginal / se.marginal))) * 2
        est.marginal <- cbind(out$marginal, se.marginal, CI.marginal, pvalue.marginal)
        colnames(est.marginal) <- c("marginal", "S.E.", "CI.lower", "CI.upper", "p.value")
      }
    }

    ## placebo test
    if (!is.null(placebo.period) & placeboTest == TRUE) {
      att.placebo <- out$att.placebo
      se.placebo <- sd(att.placebo.boot, na.rm = TRUE)
      if (quantile.CI == FALSE) {
        CI.placebo <- c(
          att.placebo - se.placebo * qnorm(1 - alpha / 2),
          att.placebo + se.placebo * qnorm(1 - alpha / 2)
        )
        CI.placebo.bound <- c(
          att.placebo - se.placebo * qnorm(1 - alpha),
          att.placebo + se.placebo * qnorm(1 - alpha)
        )
        pvalue.placebo <- (1 - pnorm(abs(att.placebo / se.placebo))) * 2
      } else {
        CI.placebo <- quantile(att.placebo.boot, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
        CI.placebo.bound <- quantile(att.placebo.boot, c(alpha, 1 - alpha), na.rm = TRUE)
        pvalue.placebo <- get.pvalue(att.placebo.boot)
      }

      est.placebo <- t(as.matrix(c(
        att.placebo,
        se.placebo,
        CI.placebo,
        pvalue.placebo,
        CI.placebo.bound
      )))
      colnames(est.placebo) <- c(
        "ATT.placebo", "S.E.",
        "CI.lower", "CI.upper", "p.value",
        "CI.lower(90%)", "CI.upper(90%)"
      )
    }

    ## carryover test
    if (!is.null(carryover.period) & carryoverTest == TRUE) {
      att.carryover <- out$att.carryover
      se.carryover <- sd(att.carryover.boot, na.rm = TRUE)
      if (quantile.CI == FALSE) {
        CI.carryover <- c(
          att.carryover - se.carryover * qnorm(1 - alpha / 2),
          att.carryover + se.carryover * qnorm(1 - alpha / 2)
        )
        CI.carryover.bound <- c(
          att.carryover - se.carryover * qnorm(1 - alpha),
          att.carryover + se.carryover * qnorm(1 - alpha)
        )
        pvalue.carryover <- (1 - pnorm(abs(att.carryover / se.carryover))) * 2
      } else {
        CI.carryover <- quantile(att.carryover.boot, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
        CI.carryover.bound <- quantile(att.carryover.boot, c(alpha, 1 - alpha), na.rm = TRUE)
        pvalue.carryover <- get.pvalue(att.carryover.boot)
      }
      est.carryover <- t(as.matrix(c(
        att.carryover, se.carryover,
        CI.carryover, pvalue.carryover,
        CI.carryover.bound
      )))
      colnames(est.carryover) <- c(
        "ATT.carryover", "S.E.",
        "CI.lower", "CI.upper", "p.value",
        "CI.lower(90%)", "CI.upper(90%)"
      )
    }

    ## group effect
    if (!is.null(group)) {
      se.group.att <- apply(group.att.boot, 1, function(vec) sd(vec, na.rm = TRUE))
      if (quantile.CI == TRUE) {
        CI.group.att <- cbind(
          c(out$group.att) - se.group.att * qnorm(1 - alpha / 2),
          c(out$group.att) + se.group.att * qnorm(1 - alpha / 2)
        )
        pvalue.group.att <- (1 - pnorm(abs(out$group.att / se.group.att))) * 2
      } else {
        CI.group.att <- t(apply(group.att.boot, 1, function(vec) quantile(vec, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)))
        pvalue.group.att <- apply(group.att.boot, 1, get.pvalue)
      }

      est.group.att <- cbind(out$group.att, se.group.att, CI.group.att, pvalue.group.att)
      colnames(est.group.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value")

      est.group.out <- list()
      for (sub.name in group.output.name) {
        subgroup.atts <- group.output.origin[[sub.name]]$att.on
        subgroup.atts.boot <- group.atts.boot[[sub.name]]
        subgroup.est.att <- NULL
        subgroup.att.bound <- NULL
        if (dim(subgroup.atts.boot)[1] > 0) {
          subgroup.se.att <- apply(subgroup.atts.boot, 1, function(vec) sd(vec, na.rm = TRUE))
          if (quantile.CI == FALSE) {
            subgroup.CI.att <- cbind(
              subgroup.atts - subgroup.se.att * qnorm(1 - alpha / 2),
              subgroup.atts + subgroup.se.att * qnorm(1 - alpha / 2)
            )
            subgroup.pvalue.att <- (1 - pnorm(abs(subgroup.atts / subgroup.se.att))) * 2
            subgroup.att.bound <- cbind(
              subgroup.atts - subgroup.se.att * qnorm(1 - alpha),
              subgroup.atts + subgroup.se.att * qnorm(1 - alpha)
            )
          } else {
            subgroup.CI.att <- t(apply(subgroup.atts.boot, 1, function(vec) quantile(vec, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)))
            subgroup.pvalue.att <- apply(subgroup.atts.boot, 1, get.pvalue)
            subgroup.att.bound <- t(apply(subgroup.atts.boot, 1, function(vec) quantile(vec, c(alpha, 1 - alpha), na.rm = TRUE)))
          }
          subgroup.est.att <- cbind(
            subgroup.atts, subgroup.se.att,
            subgroup.CI.att, subgroup.pvalue.att,
            group.output.origin[[sub.name]]$count.on
          )
          colnames(subgroup.est.att) <- c(
            "ATT", "S.E.", "CI.lower", "CI.upper",
            "p.value", "count"
          )
          rownames(subgroup.est.att) <- group.output.origin[[sub.name]]$time.on

          # for equivalence test
          colnames(subgroup.att.bound) <- c("CI.lower", "CI.upper")
          rownames(subgroup.att.bound) <- group.output.origin[[sub.name]]$time.on
        }

        subgroup.att.off.bound <- NULL
        subgroup.est.att.off <- NULL
        if (hasRevs == 1) {
          subgroup.atts.off <- group.output.origin[[sub.name]]$att.off
          subgroup.atts.off.boot <- group.atts.off.boot[[sub.name]]

          if (dim(subgroup.atts.off.boot)[1] > 0) {
            subgroup.se.att.off <- apply(subgroup.atts.off.boot, 1, function(vec) sd(vec, na.rm = TRUE))
            if (quantile.CI == FALSE) {
              subgroup.CI.att.off <- cbind(
                subgroup.atts.off - subgroup.se.att.off * qnorm(1 - alpha / 2),
                subgroup.atts.off + subgroup.se.att.off * qnorm(1 - alpha / 2)
              )
              subgroup.pvalue.att.off <- apply(subgroup.atts.off.boot, 1, get.pvalue)
              subgroup.att.off.bound <- cbind(
                subgroup.atts.off - subgroup.se.att.off * qnorm(1 - alpha),
                subgroup.atts.off + subgroup.se.att.off * qnorm(1 - alpha)
              )
            } else {
              subgroup.CI.att.off <- t(apply(subgroup.atts.off.boot, 1, function(vec) quantile(vec, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)))
              subgroup.pvalue.att.off <- apply(subgroup.atts.off.boot, 1, get.pvalue)
              subgroup.att.off.bound <- t(apply(subgroup.atts.off.boot, 1, function(vec) quantile(vec, c(alpha, 1 - alpha), na.rm = TRUE)))
            }
            subgroup.est.att.off <- cbind(
              subgroup.atts.off,
              subgroup.se.att.off,
              subgroup.CI.att.off,
              subgroup.pvalue.att.off,
              group.output.origin[[sub.name]]$count.off
            )
            colnames(subgroup.est.att.off) <- c(
              "ATT.OFF", "S.E.", "CI.lower", "CI.upper",
              "p.value", "count.off"
            )
            rownames(subgroup.est.att.off) <- group.output.origin[[sub.name]]$time.off

            colnames(subgroup.att.off.bound) <- c("CI.lower", "CI.upper")
            rownames(subgroup.att.off.bound) <- group.output.origin[[sub.name]]$time.off
          }
        }

        ## placebo test
        subgroup.est.placebo <- NULL
        if (!is.null(placebo.period) & placeboTest == TRUE) {
          subgroup.att.placebo <- group.output.origin[[sub.name]]$att.placebo
          if (length(subgroup.att.placebo) > 0) {
            subgroup.se.placebo <- sd(group.att.placebo.boot[[sub.name]], na.rm = TRUE)
            if (quantile.CI == FALSE) {
              subgroup.CI.placebo <- c(
                subgroup.att.placebo - subgroup.se.placebo * qnorm(1 - alpha / 2),
                subgroup.att.placebo + subgroup.se.placebo * qnorm(1 - alpha / 2)
              )
              subgroup.CI.placebo.bound <- c(
                subgroup.att.placebo - subgroup.se.placebo * qnorm(1 - alpha),
                subgroup.att.placebo + subgroup.se.placebo * qnorm(1 - alpha)
              )
              subgroup.pvalue.placebo <- (1 - pnorm(abs(subgroup.att.placebo / subgroup.se.placebo))) * 2
            } else {
              subgroup.CI.placebo <- quantile(group.att.placebo.boot[[sub.name]], c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
              subgroup.CI.placebo.bound <- quantile(group.att.placebo.boot[[sub.name]], c(alpha, 1 - alpha), na.rm = TRUE)
              subgroup.pvalue.placebo <- get.pvalue(group.att.placebo.boot[[sub.name]])
            }

            subgroup.est.placebo <- t(as.matrix(c(
              subgroup.att.placebo,
              subgroup.se.placebo,
              subgroup.CI.placebo,
              subgroup.pvalue.placebo,
              subgroup.CI.placebo.bound
            )))
            colnames(subgroup.est.placebo) <- c(
              "ATT.placebo", "S.E.",
              "CI.lower", "CI.upper", "p.value",
              "CI.lower(90%)", "CI.upper(90%)"
            )
          }
        }

        ## carryover test
        subgroup.est.carryover <- NULL
        if (!is.null(carryover.period) & carryoverTest == TRUE) {
          subgroup.att.carryover <- group.output.origin[[sub.name]]$att.carryover
          if (length(subgroup.att.carryover) > 0) {
            subgroup.se.carryover <- sd(group.att.carryover.boot[[sub.name]], na.rm = TRUE)
            if (quantile.CI == FALSE) {
              subgroup.CI.carryover <- c(
                subgroup.att.carryover - subgroup.se.carryover * qnorm(1 - alpha / 2),
                subgroup.att.carryover + subgroup.se.carryover * qnorm(1 - alpha / 2)
              )
              subgroup.CI.carryover.bound <- c(
                subgroup.att.carryover - subgroup.se.carryover * qnorm(1 - alpha),
                subgroup.att.carryover + subgroup.se.carryover * qnorm(1 - alpha)
              )
              subgroup.pvalue.carryover <- (1 - pnorm(abs(subgroup.att.carryover / subgroup.se.carryover))) * 2
            } else {
              subgroup.CI.carryover <- quantile(group.att.carryover.boot[[sub.name]], c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
              subgroup.CI.carryover.bound <- quantile(group.att.carryover.boot[[sub.name]], c(alpha, 1 - alpha), na.rm = TRUE)
              subgroup.pvalue.carryover <- get.pvalue(group.att.carryover.boot[[sub.name]])
            }

            subgroup.est.carryover <- t(as.matrix(c(
              subgroup.att.carryover,
              subgroup.se.carryover,
              subgroup.CI.carryover,
              subgroup.pvalue.carryover,
              subgroup.CI.carryover.bound
            )))
            colnames(subgroup.est.carryover) <- c(
              "ATT.carryover", "S.E.",
              "CI.lower", "CI.upper", "p.value",
              "CI.lower(90%)", "CI.upper(90%)"
            )
          }
        }

        est.group.out[[sub.name]] <- list(
          att.on = subgroup.est.att,
          att.on.bound = subgroup.att.bound,
          att.on.boot = group.atts.boot[[sub.name]],
          att.off = subgroup.est.att.off,
          att.off.bound = subgroup.att.off.bound,
          att.off.boot = group.atts.off.boot[[sub.name]],
          att.placebo = subgroup.est.placebo,
          att.carryover = subgroup.est.carryover
        )
      }
    }
  }

  ## storage
  result <- list(
    est.avg = est.avg,
    att.bound = att.bound,
    att.avg.boot = att.avg.boot,
    est.avg.unit = est.avg.unit,
    att.avg.unit.boot = att.avg.unit.boot,
    est.eff.calendar = est.eff.calendar,
    est.eff.calendar.fit = est.eff.calendar.fit,
    est.att = est.att,
    est.att90 = est.att90,
    att.boot = att.boot,
    att.boot.original = att.boot.original,
    att.vcov = vcov.att,
    att.count.boot = att.count.boot,
    vartype = vartype
  )
  if (keep.sims) {
    result = c(result, list(
      eff.boot = eff.boot,
      D.boot = D.boot,
      I.boot = I.boot,
      colnames.boot = colnames.boot
    ))
  }

  if (p > 0) {
    result <- c(result, list(beta.boot = beta.boot))
    result <- c(result, list(est.beta = est.beta))
    if (binary == TRUE) {
      result <- c(result, list(est.marginal = est.marginal))
    }
  }
  if (hasRevs == 1) {
    result <- c(result, list(
      est.att.off = est.att.off,
      att.off.boot = att.off.boot,
      att.off.vcov = vcov.att.off,
      att.off.bound = att.off.bound,
      att.off.count.boot = att.off.count.boot
    ))
    if (keep.sims) {
      result <- c(result, list(eff.off.boot = eff.off.boot))
    }
  }
  if (!is.null(T.on.carry)) {
    result <- c(result, list(est.carry.att = est.carry.att))
  }

  if (!is.null(balance.period)) {
    result <- c(result, list(est.balance.att = est.balance.att))
    result <- c(result, list(est.balance.avg = est.balance.avg))
    result <- c(result, list(
      balance.att.bound = balance.att.bound,
      balance.att.vcov = vcov.balance.att,
      balance.att.boot = balance.att.boot,
      balance.count.boot = balance.count.boot
    ))
    if (!is.null(placebo.period) & placeboTest == TRUE) {
      result <- c(result, list(est.balance.placebo = est.balance.placebo, balance.att.placebo.boot = balance.att.placebo.boot))
    }
  }
  if (!is.null(W)) {
    # att.avg.W.boot
    result <- c(result, list(est.avg.W = est.avg.W))
    result <- c(result, list(est.att.W = est.att.W))
    result <- c(result, list(att.W.bound = att.W.bound))
    result <- c(result, list(att.W.boot = att.on.W.boot, att.W.vcov = vcov.att.W))
    if (!is.null(placebo.period) & placeboTest == TRUE) {
      result <- c(result, list(est.placebo.W = est.placebo.W))
    }
    if (hasRevs == 1) {
      result <- c(result, list(est.att.off.W = est.att.off.W, att.off.W.bound = att.off.W.bound, att.off.W.vcov = vcov.att.off.W))
      if (!is.null(carryover.period) & carryoverTest == TRUE) {
        result <- c(result, list(est.carryover.W = est.carryover.W))
      }
    }
  }


  if (!is.null(placebo.period) & placeboTest == TRUE) {
    result <- c(result, list(est.placebo = est.placebo, att.placebo.boot = att.placebo.boot))
  }

  if (!is.null(carryover.period) & carryoverTest == TRUE) {
    result <- c(result, list(est.carryover = est.carryover, att.carryover.boot = att.carryover.boot))
  }

  if (!is.null(group)) {
    result <- c(result, list(
      est.group.att = est.group.att,
      est.group.output = est.group.out
    ))
  }




  return(c(out, result))
} ## end of boot


## jackknife se
jackknifed <- function(x, ## ols estimates
                       y,
                       alpha,
                       quantile.CI = FALSE) { ## sub-sample ols estimates)

  p <- length(x)
  N <- dim(y)[2] ## sample size

  X <- matrix(rep(c(x), N), p, N) * N
  Y <- X - y * (N - 1)

  Yvar <- apply(Y, 1, var, na.rm = TRUE)
  vn <- N - apply(is.na(y), 1, sum)

  Ysd <- sqrt(Yvar / vn) ## jackknife se

  vcov_matrix <- 1 / vn * cov(t(Y), use = "pairwise.complete.obs")

  if (quantile.CI == FALSE) {
    CI.l <- Ysd * qnorm(alpha / 2) + c(x)
    CI.u <- Ysd * qnorm(1 - alpha / 2) + c(x)
  } else {
    CI <- t(apply(y, 1, function(vec) quantile(vec, c(0.05 / 2, 1 - 0.05 / 2), na.rm = TRUE)))
    CI.l <- CI[, 1]
    CI.u <- CI[, 2]
  }


  ## wald test
  P <- NULL
  for (i in 1:p) {
    subz <- pnorm(c(x)[i] / Ysd[i])
    P <- c(P, 2 * min(1 - subz, subz))
  }

  ## P <- 2 * min(1 - pnorm(c(x)/Ysd), pnorm(c(x)/Ysd))

  out <- list(se = Ysd, CI.l = CI.l, CI.u = CI.u, P = P, vcov = vcov_matrix)

  return(out)
}
