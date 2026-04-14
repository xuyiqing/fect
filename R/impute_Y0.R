## impute_Y0() — unified Y0-imputer dispatcher for parametric bootstrap
##
## Dispatches to the correct underlying estimator based on (method, predictive).
## Pure dispatcher — introduces NO RNG calls.  try() wrapping is the caller's
## responsibility (see draw.error() and one.nonpara() in boot.R).
##
## Dispatch table:
##   gsynth  OR (ife, nevertreated)   -> fect_nevertreated(method="ife",  r=tuning)
##   ife,  notyettreated              -> fect_fe(r.cv=tuning)
##   cfe,  nevertreated               -> fect_nevertreated(method="cfe",  r=tuning)
##   cfe,  notyettreated              -> fect_cfe(r.cv=tuning)
##   mc    (any)                      -> stop()  [Phase 4 – not yet implemented]

impute_Y0 <- function(
  method,           # "gsynth", "ife", "cfe"  (dispatch key 1)
  predictive,       # "notyettreated", "nevertreated"  (dispatch key 2)

  ## Data arguments (passed through to callee)
  Y,                # TT x N outcome matrix
  X,                # TT x N x p covariate array, or NULL
  D,                # TT x N treatment indicator
  W,                # TT x N weight matrix or NULL
  I,                # TT x N observation indicator
  II,               # TT x N estimation indicator
  T.on,             # TT x N integer matrix: periods since onset

  ## Extra-FE arguments — forwarded only to cfe-family callees; ignored by fect_fe
  X.extra.FE = NULL,
  X.Z        = NULL,
  X.Q        = NULL,
  X.gamma    = NULL,
  X.kappa    = NULL,
  Zgamma.id  = NULL,
  kappaQ.id  = NULL,

  ## Tuning and fit control
  tuning,           # scalar: r.cv (pre-selected rank). Forwarded as r= or r.cv= per callee.
  boot = 1,         # always 1 for bootstrap refits; kept for signature symmetry

  ## Aggregation sequences (forwarded to callee)
  T.on.balance     = NULL,
  balance.period   = NULL,
  placeboTest      = 0,
  placebo.period   = NULL,
  carryoverTest    = 0,
  carryover.period = NULL,
  calendar.enp.seq = NULL,
  time.on.seq      = NULL,
  time.off.seq     = NULL,
  time.on.seq.W    = NULL,
  time.off.seq.W   = NULL,
  time.on.seq.group = NULL,
  time.off.seq.group = NULL,
  time.on.balance.seq = NULL,

  ## Estimation control (forwarded to callee)
  force,
  hasRevs      = 1,
  tol,
  max.iteration = 1000,
  norm.para    = NULL,
  group.level  = NULL,
  group        = NULL
) {

  if (method == "gsynth" || (method == "ife" && predictive == "nevertreated")) {
    ## Branch 1: gsynth and ife+nevertreated both call fect_nevertreated(method="ife")
    fect_nevertreated(
      Y              = Y,
      X              = X,
      D              = D,
      W              = W,
      I              = I,
      II             = II,
      T.on           = T.on,
      T.on.balance   = T.on.balance,
      balance.period = balance.period,
      hasRevs        = hasRevs,
      force          = force,
      r              = tuning,
      CV             = 0,
      boot           = boot,
      placeboTest    = placeboTest,
      placebo.period = placebo.period,
      carryover.period = carryover.period,
      carryoverTest  = carryoverTest,
      calendar.enp.seq = calendar.enp.seq,
      time.on.seq    = time.on.seq,
      time.off.seq   = time.off.seq,
      time.on.seq.W  = time.on.seq.W,
      time.off.seq.W = time.off.seq.W,
      time.on.seq.group  = time.on.seq.group,
      time.off.seq.group = time.off.seq.group,
      time.on.balance.seq = time.on.balance.seq,
      norm.para      = norm.para,
      tol            = tol,
      max.iteration  = max.iteration,
      group.level    = group.level,
      group          = group,
      method         = "ife"
    )

  } else if (method == "ife" && predictive == "notyettreated") {
    ## Branch 2: ife+notyettreated -> fect_fe
    fect_fe(
      Y              = Y,
      X              = X,
      D              = D,
      W              = W,
      I              = I,
      II             = II,
      T.on           = T.on,
      T.on.balance   = T.on.balance,
      balance.period = balance.period,
      hasRevs        = hasRevs,
      force          = force,
      r.cv           = tuning,
      boot           = boot,
      placeboTest    = placeboTest,
      placebo.period = placebo.period,
      carryover.period = carryover.period,
      carryoverTest  = carryoverTest,
      calendar.enp.seq = calendar.enp.seq,
      time.on.seq    = time.on.seq,
      time.off.seq   = time.off.seq,
      time.on.seq.W  = time.on.seq.W,
      time.off.seq.W = time.off.seq.W,
      time.on.seq.group  = time.on.seq.group,
      time.off.seq.group = time.off.seq.group,
      time.on.balance.seq = time.on.balance.seq,
      tol            = tol,
      max.iteration  = max.iteration,
      norm.para      = norm.para,
      group.level    = group.level,
      group          = group
    )

  } else if (method == "cfe" && predictive == "nevertreated") {
    ## Branch 3: cfe+nevertreated -> fect_nevertreated(method="cfe")
    fect_nevertreated(
      Y              = Y,
      X              = X,
      D              = D,
      W              = W,
      I              = I,
      II             = II,
      T.on           = T.on,
      T.on.balance   = T.on.balance,
      balance.period = balance.period,
      hasRevs        = hasRevs,
      force          = force,
      r              = tuning,
      CV             = 0,
      boot           = boot,
      placeboTest    = placeboTest,
      placebo.period = placebo.period,
      carryover.period = carryover.period,
      carryoverTest  = carryoverTest,
      calendar.enp.seq = calendar.enp.seq,
      time.on.seq    = time.on.seq,
      time.off.seq   = time.off.seq,
      time.on.seq.W  = time.on.seq.W,
      time.off.seq.W = time.off.seq.W,
      time.on.seq.group  = time.on.seq.group,
      time.off.seq.group = time.off.seq.group,
      time.on.balance.seq = time.on.balance.seq,
      norm.para      = norm.para,
      tol            = tol,
      max.iteration  = max.iteration,
      group.level    = group.level,
      group          = group,
      method         = "cfe",
      X.extra.FE    = X.extra.FE,
      X.Z           = X.Z,
      X.Q           = X.Q,
      X.gamma       = X.gamma,
      X.kappa       = X.kappa,
      Zgamma.id     = Zgamma.id,
      kappaQ.id     = kappaQ.id
    )

  } else if (method == "cfe" && predictive == "notyettreated") {
    ## Branch 4: cfe+notyettreated -> fect_cfe
    fect_cfe(
      Y              = Y,
      X              = X,
      D              = D,
      W              = W,
      X.extra.FE    = X.extra.FE,
      X.Z           = X.Z,
      X.Q           = X.Q,
      X.gamma       = X.gamma,
      X.kappa       = X.kappa,
      Zgamma.id     = Zgamma.id,
      kappaQ.id     = kappaQ.id,
      I              = I,
      II             = II,
      T.on           = T.on,
      T.on.balance   = T.on.balance,
      balance.period = balance.period,
      hasRevs        = hasRevs,
      force          = force,
      r.cv           = tuning,
      boot           = boot,
      placeboTest    = placeboTest,
      placebo.period = placebo.period,
      carryover.period = carryover.period,
      carryoverTest  = carryoverTest,
      calendar.enp.seq = calendar.enp.seq,
      time.on.seq    = time.on.seq,
      time.off.seq   = time.off.seq,
      time.on.seq.W  = time.on.seq.W,
      time.off.seq.W = time.off.seq.W,
      time.on.seq.group  = time.on.seq.group,
      time.off.seq.group = time.off.seq.group,
      time.on.balance.seq = time.on.balance.seq,
      tol            = tol,
      max.iteration  = max.iteration,
      norm.para      = norm.para,
      group.level    = group.level,
      group          = group
    )

  } else {
    stop(paste0(
      "impute_Y0: unsupported combination method='", method,
      "', predictive='", predictive, "'. ",
      "MC support is planned for Phase 4."
    ))
  }
}

## partition_controls() — internal helper for K=2 cross-fitting debiasing
##
## Partitions a control-unit index vector into K=2 disjoint halves.
## Half-A (smaller or equal): used as the factor-estimation pool in Loop 1.
## Half-B (larger or equal): used as the OOS-residual pool in Loop 2.
##
## Uses the *current* RNG state — no internal set.seed() call.
## When split_residuals=FALSE this function is never called, so the FALSE path
## advances the RNG by exactly zero extra calls (parity contract preserved).
##
## Args:
##   id.co  — integer vector of control column indices, length Nco (>= 4 required)
##   K      — integer scalar, must equal 2L for this release
##
## Returns: list with elements $A (half-A indices) and $B (half-B indices),
##          each sorted to preserve column order.

partition_controls <- function(id.co, K = 2L) {
  Nco <- length(id.co)

  ## --- Input validation ---
  if (Nco < 4L) {
    stop(
      "split_residuals requires at least 4 control units (Nco = ", Nco,
      "). Reduce K or disable split_residuals."
    )
  }
  if (!identical(K, 2L)) {
    stop("partition_controls: only K=2 is supported in this release.")
  }

  ## --- Partition sizes ---
  Nco_A <- floor(Nco / 2L)
  ## Nco_B = Nco - Nco_A  (= ceiling(Nco/2), the larger half)

  ## --- Random assignment (uses current RNG state) ---
  idx   <- sample(Nco)                          # random permutation of 1:Nco
  A_idx <- idx[seq_len(Nco_A)]
  B_idx <- idx[seq.int(Nco_A + 1L, Nco)]

  ## --- Sort within each half to preserve column order ---
  id.co_A <- id.co[sort(A_idx)]
  id.co_B <- id.co[sort(B_idx)]

  list(A = id.co_A, B = id.co_B)
}
