#' Over-identification test for causal moderation (cm)
#'
#' Implements the over-identification test described in the slide:
#' run two auxiliary regressions of residualized outcomes on the moderator
#' and covariates (including quadratic and interaction terms) with unit and
#' time fixed effects, then test using \(n \times R^2 \sim \chi^2(df)\),
#' where df is the number of *nonlinear* (quadratic + interaction) terms that
#' remain in the fitted model.
#'
#' @param x A `fect` object. Must be estimated with `cm=TRUE` so that `x$est.cm` exists.
#' @param moderator Character scalar. Name of the moderator \(M_{it}\) (must be in `x$X`).
#' @param covariates Optional character vector. Names of other covariates \(X_{it}\).
#'   Default is all covariates in `x$X` except `moderator`.
#' @param quadratic Logical. Include quadratic terms \(M_{it}^2\) and \(X_{it}^2\). Default TRUE.
#' @param interaction Logical. Include interaction terms \(M_{it} \times X_{it}\). Default TRUE.
#'
#' @return A list with elements `e1` and `e0`, each containing:
#'   - `n`: sample size used in the auxiliary regression
#'   - `r2`: R-squared of the auxiliary regression
#'   - `stat`: \(n \times R^2\)
#'   - `df`: degrees of freedom (# nonlinear terms kept)
#'   - `p`: p-value from chi-square test
#'   - `model`: the fitted `fixest` model
#'
#' @export
fect_iden <- function(x,
                      moderator,
                      covariates = NULL,
                      quadratic = TRUE,
                      interaction = TRUE) {
  if (is.null(x$est) || is.null(x$est$fit)) {
    stop("`x$est$fit` is required to compute g0.\n")
  }
  if (is.null(x$est.cm) || is.null(x$est.cm$fit)) {
    stop("This test requires `cm=TRUE` so that `x$est.cm$fit` exists.\n")
  }
  if (is.null(x$Y.dat) || is.null(x$D.dat)) {
    stop("`x$Y.dat` and `x$D.dat` are required.\n")
  }
  if (is.null(x$id) || is.null(x$rawtime)) {
    stop("`x$id` and `x$rawtime` are required.\n")
  }
  if (is.null(x$X) || length(x$X) == 0) {
    stop("`x$X` (covariate names) is required.\n")
  }

  if (!is.character(moderator) || length(moderator) != 1) {
    stop("`moderator` must be a single character string.\n")
  }
  if (!moderator %in% x$X) {
    stop("`moderator` must be one of `x$X`.\n")
  }

  if (is.null(covariates)) {
    covariates <- setdiff(x$X, moderator)
  }
  if (!is.character(covariates)) {
    stop("`covariates` must be a character vector.\n")
  }
  if (length(covariates) > 0 && any(!covariates %in% x$X)) {
    stop("All `covariates` must be in `x$X`.\n")
  }

  # Extract covariate array (T x N x p) from the fect object.
  # In `fect` objects there are typically two entries named "X":
  # - the first is `x$X` (character vector of covariate names)
  # - the second is the actual covariate array (T x N x p)
  X.arr <- NULL
  X.cands <- x[names(x) == "X"]
  if (length(X.cands) >= 2) {
    # Using `[2]$X` (not `[[2]]$X`) works whether the second element is a list or an array.
    X.arr <- X.cands[2]$X
  } else if (is.array(x$X)) {
    # Fallback: some internal objects may store the array directly at x$X
    X.arr <- x$X
  }
  if (is.null(X.arr) || length(dim(X.arr)) != 3) {
    stop("Cannot locate covariate array in `x`.\n")
  }

  TT <- nrow(x$Y.dat)
  N <- ncol(x$Y.dat)
  if (nrow(x$D.dat) != TT || ncol(x$D.dat) != N) {
    stop("Dimensions of `x$Y.dat` and `x$D.dat` do not match.\n")
  }

  g0 <- x$est$fit
  g1 <- x$est.cm$fit
  if (!is.matrix(g0) || !is.matrix(g1) || any(dim(g0) != c(TT, N)) || any(dim(g1) != c(TT, N))) {
    stop("`x$est$fit` and `x$est.cm$fit` must be T-by-N matrices matching `x$Y.dat`.\n")
  }

  Y <- x$Y.dat
  D <- x$D.dat
  I <- if (!is.null(x$I.dat)) x$I.dat else matrix(1, TT, N)

  # Residual-like objects defined in the slide:
  # e1_it = D_it (Y_it - g1_it), on treated cells
  # e0_it = (1-D_it) (Y_it - g0_it), on control cells
  e1 <- D * (Y - g1)
  e0 <- (1 - D) * (Y - g0)

  # Long-form indexing consistent with c(matrix): time varies fastest within unit.
  unit <- rep(x$id, each = TT)
  time <- rep(x$rawtime, times = N)

  # Extract moderator and covariates into vectors aligned with c(Y).
  m_idx <- which(x$X == moderator)[1]
  M.vec <- c(X.arr[, , m_idx])

  X.dat <- list()
  if (length(covariates) > 0) {
    for (nm in covariates) {
      idx <- which(x$X == nm)[1]
      X.dat[[nm]] <- c(X.arr[, , idx])
    }
  }

  Y.vec <- c(Y)
  D.vec <- c(D)
  I.vec <- c(I)
  g0.vec <- c(g0)
  g1.vec <- c(g1)
  e1.vec <- c(e1)
  e0.vec <- c(e0)

  base_df <- data.frame(
    unit = as.factor(unit),
    time = as.factor(time),
    M = M.vec
  )
  names(base_df)[names(base_df) == "M"] <- moderator
  if (length(X.dat) > 0) {
    for (nm in names(X.dat)) {
      base_df[[nm]] <- X.dat[[nm]]
    }
  }

  build_formula <- function(lhs) {
    linear_terms <- c(moderator, covariates)
    nonlinear_terms <- c()

    if (isTRUE(quadratic)) {
      nonlinear_terms <- c(nonlinear_terms, paste0("I(", moderator, "^2)"))
      if (length(covariates) > 0) {
        nonlinear_terms <- c(nonlinear_terms, paste0("I(", covariates, "^2)"))
      }
    }
    if (isTRUE(interaction) && length(covariates) > 0) {
      nonlinear_terms <- c(nonlinear_terms, paste0(moderator, ":", covariates))
    }

    rhs_terms <- c(linear_terms, nonlinear_terms)
    rhs_terms <- rhs_terms[nzchar(rhs_terms)]
    rhs <- if (length(rhs_terms) > 0) paste(rhs_terms, collapse = " + ") else "1"
    stats::as.formula(paste0(lhs, " ~ ", rhs, " | unit + time"))
  }

  r2_of <- function(mod) {
    # Prefer fixest fitstat if available; fall back to r2().
    r2v <- tryCatch({
      fs <- fixest::fitstat(mod, "r2")
      as.numeric(fs$r2)
    }, error = function(e) NA_real_)
    if (is.na(r2v)) {
      r2v <- tryCatch(as.numeric(fixest::r2(mod, type = "r2")), error = function(e) NA_real_)
    }
    r2v
  }

  df_nonlinear <- function(mod) {
    cn <- names(stats::coef(mod))
    if (is.null(cn) || length(cn) == 0) return(0L)

    keep <- rep(FALSE, length(cn))
    # quadratic terms
    if (isTRUE(quadratic)) {
      keep <- keep | grepl(paste0("^I\\(", moderator, "\\^2\\)$"), cn)
      if (length(covariates) > 0) {
        for (nm in covariates) {
          keep <- keep | grepl(paste0("^I\\(", nm, "\\^2\\)$"), cn)
        }
      }
    }
    # interaction terms
    if (isTRUE(interaction) && length(covariates) > 0) {
      for (nm in covariates) {
        keep <- keep | (cn == paste0(moderator, ":", nm)) | (cn == paste0(nm, ":", moderator))
      }
    }
    sum(keep)
  }

  run_one <- function(resp_vec, keep_idx, lhs_name) {
    dat <- base_df[keep_idx, , drop = FALSE]
    dat[[lhs_name]] <- resp_vec[keep_idx]

    fml <- build_formula(lhs_name)
    mod <- fixest::feols(fml, data = dat, warn = FALSE)

    n <- stats::nobs(mod)
    r2 <- r2_of(mod)
    df <- df_nonlinear(mod)
    stat <- n * r2
    p <- if (is.finite(stat) && df > 0) (1 - stats::pchisq(stat, df = df)) else NA_real_

    list(n = n, r2 = r2, stat = stat, df = df, p = p, model = mod)
  }

  # Restrict to observed cells (I==1) and appropriate treatment status.
  keep_common <- which(I.vec == 1 & !is.na(Y.vec))
  keep_e1 <- intersect(keep_common, which(D.vec == 1 & !is.na(g1.vec) & !is.na(e1.vec)))
  keep_e0 <- intersect(keep_common, which(D.vec == 0 & !is.na(g0.vec) & !is.na(e0.vec)))

  out <- list(
    e1 = run_one(e1.vec, keep_e1, "e1"),
    e0 = run_one(e0.vec, keep_e0, "e0"),
    call = match.call()
  )
  class(out) <- "fect_iden"
  out
}


