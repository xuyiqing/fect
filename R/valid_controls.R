## valid_controls() — control screening helper for parametric bootstrap
##
## Returns the subset of id.co indices that pass the pre/post-period
## observation threshold.  Used by fect_boot() before draw.error() runs.
##
## Arguments:
##   out        – original-fit result list (contains $D, $I, $r.cv)
##   method     – character scalar: "gsynth", "ife", "cfe"  (mc is blocked upstream)
##   predictive – character scalar: "notyettreated" or "nevertreated"
##   force      – integer scalar: 0, 1, 2, or 3

valid_controls <- function(out, method, predictive, force) {

  if (method == "mc") {
    stop("valid_controls: MC not supported in this release.")
  }

  id.co  <- which(colSums(out$D) == 0)
  id.tr  <- which(colSums(out$D) > 0)
  I.co   <- out$I[, id.co, drop = FALSE]
  TT     <- nrow(out$D)

  T0.ub      <- colSums(out$D[, id.tr, drop = FALSE] == 0)
  T0.ub.min  <- min(T0.ub)
  max.T0.ub  <- max(T0.ub)

  co.pre  <- colSums(I.co[seq_len(T0.ub.min), , drop = FALSE])
  co.post <- colSums(I.co[seq.int(max.T0.ub + 1L, TT), , drop = FALSE])

  threshold <- out$r.cv + as.integer(force %in% c(1, 3))

  valid.co <- id.co[(co.pre >= threshold) & (co.post >= 1L)]
  return(valid.co)
}
