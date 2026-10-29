##########
## Print
##########
## a fect object

print.fect <- function(x,
                       switch.on = TRUE,
                       switch.off = FALSE,
                       time.on.lim = NULL,
                       time.off.lim = NULL,
                       ...) {


    if (!is.null(x$effect.est.avg)) {
      cat("\nOverall cumulative effect:\n")
      print(x$effect.est.avg, digits = 4)
      if (!is.null(x$effect.est.att)) {
        cat("\nPeriod-by-period cumulative effect:\n")
        print(x$effect.est.att, digits = 4)
      }
      return()
    }
    cat("Call:\n")
    print(x$call, digits = 4)

    ## Estimator + Fixed effects + Cluster SE lines (D8 in design memo).
    ## Surfaces the actual FE composition so users can verify their model
    ## without re-reading the call args. force is stored as integer 0/1/2/3.
    .fe_label <- function(force_int, index, group.fe) {
        parts <- character(0)
        if (length(index) >= 1 && force_int %in% c(1, 3)) {
            parts <- c(parts, paste0(index[1], " (unit)"))
        }
        if (length(index) >= 2 && force_int %in% c(2, 3)) {
            parts <- c(parts, paste0(index[2], " (time)"))
        }
        if (!is.null(group.fe) && length(group.fe) > 0) {
            parts <- c(parts, paste(group.fe, collapse = " + "))
        }
        if (length(parts) == 0) "none" else paste(parts, collapse = " + ")
    }
    fe.line <- .fe_label(x$force,
                         if (!is.null(x$call$index)) eval(x$call$index) else NULL,
                         x$group.fe)
    cat("\nEstimator:    ", x$method, "\n", sep = "")
    cat("Fixed effects: ", fe.line, "\n", sep = "")
    ## The parametric bootstrap and the jackknife ignore cl, so the SEs of
    ## such fits are not clustered; say so only when they are.
    if (!is.null(x$cl.label) &&
        !isTRUE(x$vartype %in% c("parametric", "jackknife"))) {
        cat("Cluster SE:   ", x$cl.label, "\n", sep = "")
    }
    .print_convergence(x$convergence)

    if (switch.on == TRUE) {
        if (!is.null(time.on.lim)) {

            if (is.numeric(time.on.lim)==FALSE) {
                stop("Some element in \"time.on.lim\" is not numeric.")
            } else {
                if (length(time.on.lim)!=2) {
                    stop("time.on.lim must be of length 2.")
                }
            }

            seq.on.min <- min(which(x$time >= time.on.lim[1]))
            seq.on.max <- max(which(x$time <= time.on.lim[2]))
            seq.on <- seq.on.min:seq.on.max
        } else {
            seq.on <- 1:length(x$time)
        }
    }

    if (switch.off == TRUE & is.null(x$att.off) == FALSE) {
        if (!is.null(time.off.lim)) {

            if (is.numeric(time.off.lim)==FALSE) {
                stop("Some element in \"time.off.lim\" is not numeric.")
            } else {
                if (length(time.off.lim)!=2) {
                    stop("time.off.lim must be of length 2.")
                }
            }

            seq.off.min <- min(which(x$time.off >= time.off.lim[1]))
            seq.off.max <- max(which(x$time.off <= time.off.lim[2]))
            seq.off <- seq.off.min:seq.off.max
        } else {
            seq.off <- 1:length(x$time.off)
        }
    }

    if (is.null(x$est.avg) == TRUE) { # no uncertainties
        cat("\nATT:\n")
        att.out <- rbind.data.frame(x$att.avg, x$att.avg.unit)
        colnames(att.out) <- c("ATT")
        rownames(att.out) <- c(
            "Tr obs. equally weighted",
            "Tr units equally weighted")
        print(att.out, digits = 4)
        # if (switch.on == TRUE) {
        #     cat("\n   ~ by Period:\n")
        #     print(x$att[seq.on], digits = 4)
        # }
        # if (switch.off == TRUE & is.null(x$att.off) == FALSE) {
        #     cat("\n   ~ Switch-off by Period:\n")
        #     print(x$att.off[seq.off], digits = 4)
        # }
        if (is.null(x$X) == FALSE) {
            cat("\nCovariates:\n")
            print(x$beta, digits = 4)
        }
        cat("\nUncertainty estimates not available.\n")
    } else {
        cat("\nATT:\n")
        att.out <- rbind.data.frame(c(x$est.avg), c(x$est.avg.unit))
        colnames(att.out) <- c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value")
        ## When W is supplied AND enters the aggregation, the obs-level row
        ## reports the sample-weighted aggregate. Otherwise (no W, or W only
        ## entered the outcome-model fit), the row is the unweighted average.
        first.row <- if (isTRUE(x$W.in.agg)) {
            "Tr obs sample-weighted (W)"
        } else {
            "Tr obs equally weighted"
        }
        rownames(att.out) <- c(
            first.row,
            "Tr units equally weighted")
        print(att.out, digits = 4)
        # if (switch.on == TRUE) {
        #     cat("\n   ~ Switch-on by Period:\n")
        #     print(x$est.att[seq.on,], digits = 4)
        # }
        # if (switch.off == TRUE & is.null(x$att.off) == FALSE) {
        #     cat("\n   ~ Switch-on by Period:\n")
        #     print(x$est.att.off[seq.off,], digits = 4)
        # }
        if (is.null(x$X) == FALSE) {
            cat("\nCovariates:\n")
            print(x$est.beta, digits = 4)
        }
    }

    if (!is.null(x$est.placebo)) {
        cat("\nPlacebo effect for pre-treatment periods:\n")
        print(x$est.placebo, digits = 4)
    }
    ## v2.4.8: flag a bootstrap whose draws break down (see ?fect, boot.diag)
    if (isTRUE(x$boot.diag$breakdown)) {
        .print_boot_diag(x$boot.diag)
    }
}

## One line on how the IFE-EM main fit stopped (x$convergence, v2.4.8).
.print_convergence <- function(conv) {
    if (is.null(conv) || is.null(conv$niter)) return(invisible(NULL))
    rc <- conv$rel.change
    rc.txt <- if (is.null(rc) || is.na(rc)) "not available" else formatC(rc, format = "e", digits = 3)
    cap <- if (is.null(conv$iter.cap)) conv$max.iteration else conv$iter.cap
    if (isTRUE(conv$converged)) {
        cat("EM convergence: stopped by tol after ", conv$niter, " rounds ",
            "(final relative change ", rc.txt, ", tol = ", format(conv$tol), ")\n", sep = "")
    } else if (isTRUE(cap == conv$max.iteration)) {
        cat("EM convergence: NOT converged; stopped at max.iteration = ", cap,
            " (final relative change ", rc.txt, ", tol = ", format(conv$tol), ")\n", sep = "")
    } else {
        ## fe_ad_iter() (no factors, no covariates) caps itself at 500 rounds
        cat("EM convergence: NOT converged; stopped at the routine's own cap of ", cap,
            " rounds, not max.iteration = ", conv$max.iteration,
            " (final relative change ", rc.txt, ", tol = ", format(conv$tol), ")\n", sep = "")
    }
    invisible(NULL)
}

## The bootstrap breakdown check on the average effect (x$boot.diag, v2.4.8).
.print_boot_diag <- function(bd, always = FALSE) {
    if (is.null(bd)) return(invisible(NULL))
    cat("\nBootstrap draws of the average effect (", bd$n.used, " of ", bd$nboots, " used):\n", sep = "")
    cat("  SD of draws (the S.E. above): ", format(signif(bd$sd, 4)),
        ";  IQR/1.349: ", format(signif(bd$iqr.sd, 4)),
        ";  ratio: ", if (is.na(bd$sd.ratio)) "not defined" else format(signif(bd$sd.ratio, 3)), "\n", sep = "")
    cat("  draws more than 20 IQRs from the median: ",
        if (is.na(bd$n.extreme)) "not defined (IQR = 0)" else bd$n.extreme,
        ";  failed draws dropped: ", bd$n.failed, "\n", sep = "")
    if (!is.null(bd$n.sim)) {
        ## the first-step pool; with para.error = "ar" a huge draw here widens
        ## every replicate while the ratio above stays near 1
        cat("  simulated-error draws (first step): ", bd$n.sim, ";  failed: ", bd$n.sim.failed,
            ";  more than 20 IQRs from the median: ",
            if (is.na(bd$n.sim.extreme)) "not defined (IQR = 0)" else bd$n.sim.extreme,
            "\n", sep = "")
    }
    if (isTRUE(bd$n.used < bd$min.draws)) {
        cat("  Fewer than ", bd$min.draws, " usable draws: the breakdown check is not applied.\n", sep = "")
    }
    if (isTRUE(bd$breakdown)) {
        cat("  The SD is more than ", bd$threshold, " times IQR/1.349: the bootstrap has broken down ",
            "and the S.E. is not reliable.\n", sep = "")
    }
    invisible(NULL)
}

## summary() of a fect fit (v2.4.8): what print() shows, plus the bootstrap
## check on the average effect whenever it was computed.
summary.fect <- function(object, ...) {
    structure(list(object = object), class = "summary.fect")
}

print.summary.fect <- function(x, ...) {
    fit <- x$object
    print.fect(fit, ...)
    if (!is.null(fit$boot.diag) && !isTRUE(fit$boot.diag$breakdown)) {
        .print_boot_diag(fit$boot.diag)
    }
    invisible(x)
}
