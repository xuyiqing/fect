# Code adapted from the 'HonestDiD' package, version 0.2.8, by Ashesh Rambachan and Jonathan Roth,
# through 'HonestDiDFEct' 1.0.0 by Ziyi Liu (github.com/lzy318/HonestDiDFEct, commit 29c831d).
# Original file: R/fect_entry.R
# Copyright (c) 2023 HonestDiD authors. MIT licence; the full notice is in inst/COPYRIGHTS.
# Copied into fect 2.5.1 without change to any numerical line; the edits are listed in the fect repository history.
#
# Entry step for fect's placebo output.
#
# fect reports coefficients for consecutive placebo periods ending at s = 0, then
# for the post periods s = 1, 2, ..., with no reference period, and numPrePeriods
# counts s = 0. HonestDiD's engine expects coefficients relative to a reference
# period that is left out of the vector. .fectEntry() turns the first form into the
# second. It differences from s = 0, which removes that coordinate, and then moves
# the post block so that the estimate of theta is zero. The move changes theta by
# $shift and leaves the bias path alone, so the accepted set moves by exactly $shift
# and .fectAddShift() puts it back. Without the move, the engine's default grids
# would be centred on zero rather than on the estimate.

.stopIfNotFectPath <- function(betahat, numPrePeriods, numPostPeriods, l_vec) {
  if (numPrePeriods < 2 || numPostPeriods < 1) {
    base::stop(base::sprintf(paste0("fect_sens() needs at least two placebo periods ending at s = 0 and one post ",
                                    "period, but got numPrePeriods = %d and numPostPeriods = %d. Use, for example, ",
                                    "fect(placeboTest = TRUE, placebo.period = c(-2, 0))."),
                             numPrePeriods, numPostPeriods), call. = FALSE)
  }

  nm <- base::names(betahat)
  if (base::is.null(nm) && base::NCOL(betahat) == 1) nm <- base::rownames(betahat)
  if (!base::is.null(nm)) {
    s    <- base::suppressWarnings(base::as.numeric(nm))
    pre  <- s[base::seq_len(numPrePeriods)]
    post <- s[numPrePeriods + base::seq_len(numPostPeriods)]
    if (base::anyNA(s) || pre[numPrePeriods] != 0 || base::any(base::diff(pre) != 1) ||
        base::any(post != base::seq_len(numPostPeriods))) {
      base::stop(base::sprintf(paste0("betahat must hold consecutive placebo periods ending at s = 0, then post ",
                                      "periods 1, 2, ..., but its names are: %s. Use a fect() run with ",
                                      "placeboTest = TRUE and placebo.period = c(-2, 0) or similar, and set ",
                                      "post.periods in fect_sens(). A run with loo = TRUE does not work here."),
                               base::paste(nm, collapse = ", ")), call. = FALSE)
    }
  }

  l <- base::as.numeric(l_vec)
  if (base::all(l == 0)) {
    base::stop("l_vec is all zeros, so it does not define a target.", call. = FALSE)
  }
  if (numPostPeriods == 1 && l != 1) {
    base::stop("With one post period, l_vec must equal 1.", call. = FALSE)
  }
}

.fectEntry <- function(betahat, sigma, numPrePeriods, numPostPeriods, l_vec) {
  .stopIfNotConformable(betahat, sigma, numPrePeriods, numPostPeriods, l_vec)
  .stopIfNotFectPath(betahat, numPrePeriods, numPostPeriods, l_vec)

  size <- numPrePeriods + numPostPeriods
  J <- base::diag(size)[-numPrePeriods, , drop = FALSE]
  J[, numPrePeriods] <- J[, numPrePeriods] - 1

  bd <- base::as.numeric(J %*% base::as.numeric(betahat))
  Sd <- J %*% base::unname(sigma) %*% base::t(J)

  l     <- base::as.numeric(l_vec)
  post  <- numPrePeriods:(size - 1)
  shift <- base::sum(l * bd[post])
  bd[post] <- bd[post] - shift * l / base::sum(l^2)

  base::list(betahat = bd, sigma = (Sd + base::t(Sd)) / 2,
             numPrePeriods = numPrePeriods - 1L, shift = shift)
}

.fectAddShift <- function(Results, shift) {
  Results$lb <- Results$lb + shift
  Results$ub <- Results$ub + shift
  base::return(Results)
}
