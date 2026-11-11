## Coverage and width of full conformal vs Jackknife+ (DID and SC rules),
## one treated unit, one target date. usage: Rscript sim.R <reps> <cores>
suppressMessages(library(fect))
args <- commandArgs(TRUE); R <- as.integer(args[1]); cores <- as.integer(args[2])
T0 <- 20; gamma <- 0.1; tau <- 2
draw <- function(J, design) {
  N <- J + 1; TT <- T0 + 1
  if (design == "additive") {
    Y <- outer(cumsum(rnorm(TT, 0, .5)), rep(1, N)) + matrix(rnorm(N), TT, N, byrow = TRUE) + matrix(rnorm(TT * N), TT, N)
  } else {
    F <- cbind(cumsum(rnorm(TT, 0, .3)), sin(2 * pi * (1:TT) / 8))
    L <- matrix(runif(N * 2, 0, 1), N, 2)
    if (design == "factor_out") L[1, ] <- c(1.6, 1.6)   # treated outside the controls' range
    Y <- F %*% t(L) * 3 + matrix(rnorm(TT * N, 0, .5), TT, N)
  }
  Y[TT, 1] <- Y[TT, 1] + tau
  list(Ypre = Y[1:T0, , drop = FALSE], ytg = Y[TT, ])
}
one <- function(r) {
  set.seed(1e4 + r); res <- list()
  for (J in c(19, 39)) for (design in c("additive", "factor_in", "factor_out")) {
    dd <- draw(J, design)
    for (rn in c("did", "sc")) {
      fun <- fect:::.fc_resolve_rule(rn)$fun
      cf <- fect:::.fc_coefs(dd$Ypre, dd$ytg, 1L, fun)
      v <- fect:::.fc_invert_affine(cf$a, cf$b, 1L, gamma)
      cov.fc <- any(v$components$lower <= tau & tau <= v$components$upper)
      jp <- fect:::.fc_jplus(dd$Ypre, dd$ytg, 1L, fun, gamma)
      res[[length(res) + 1]] <- data.frame(rep = r, J = J, design = design, rule = rn,
        fc.cover = cov.fc, fc.width = v$CI.upper - v$CI.lower, fc.disc = v$disconnected,
        jp.cover = jp$CI.lower <= tau & tau <= jp$CI.upper, jp.width = jp$CI.upper - jp$CI.lower)
    }
  }
  do.call(rbind, res)
}
x <- do.call(rbind, parallel::mclapply(seq_len(R), one, mc.cores = cores))
saveRDS(x, "out/sim.rds")
a <- aggregate(cbind(fc.cover, jp.cover, fc.disc) ~ J + design + rule, x, mean)
w <- aggregate(cbind(fc.width, jp.width) ~ J + design + rule, x, median)
print(merge(a, w), digits = 3)
cat("DONE\n")
