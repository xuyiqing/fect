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
        rownames(att.out) <- c(
            "Tr obs equally weighted",
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
}
