#include "fect.h"

// sub-functions for complex fe

/* Obtain additive fe for ub data; assume r>0 p>0*/
// [[Rcpp::export]]
List cfe_iter(const arma::cube& XX, const arma::mat& xxinv, const arma::cube& X_extra_FE,
              const arma::cube& X_Z, const arma::cube& X_Q, const arma::cube& X_gamma, const arma::cube& X_kappa,
              Rcpp::List Zgamma_id, Rcpp::List kappaQ_id,
              const arma::mat& Y,
              const arma::mat& Y0, const arma::mat& I, const arma::mat& W,
              const arma::mat& beta0, int force,
              int r, double tolerate,
              int max_iter = 1000) {
  int T = Y.n_rows;
  int N = Y.n_cols;
  int p = XX.n_slices;
  double dif = 1.0;
  int niter = 0;
  int validF = 1;
  int use_weight;
  int r_burnin;
  int d;
  if (T <= N) {
    d = T;
  } else {
    d = N;
  }

  if (Y.n_rows == W.n_rows && Y.n_cols == W.n_cols) {
    use_weight = 1;
  } else {
    use_weight = 0;
  }

  arma::mat beta(p, 1, arma::fill::zeros);
  if (beta0.n_rows == p) {
    beta = beta0;
  }
  double mu = 0;
  arma::mat VNT(r, r);
  arma::mat FE_inter_use(T, N, arma::fill::zeros); // ife
  arma::mat fit(T, N, arma::fill::zeros);
  arma::mat fit_old(T, N, arma::fill::ones);
  arma::mat U(T, N, arma::fill::zeros);
  arma::mat alpha(N, 1, arma::fill::zeros);
  arma::mat xi(T, 1, arma::fill::zeros);
  arma::mat YY = Y;
  arma::mat e(T, N, arma::fill::zeros); // residual

  List ife_inner;
  List pf;
  arma::mat F;
  arma::mat L;

  fit = Y0;
  fit_old = fit;

  arma::mat covar_fit(T, N, arma::fill::zeros);
  for (int i = 0; i < p; i++) {
    covar_fit = covar_fit + XX.slice(i) * beta(i);
  }
  arma::mat FE = fit - covar_fit;

  // Extra FE part
  arma::mat gamma;
  arma::mat kappa;
  arma::mat extra_FE_fit(T, N, arma::fill::zeros);

  int stop_burnin = 0;
  while (dif > tolerate && niter <= max_iter) {

    YY = E_adj(Y, fit, I); // e-step: expectation

    // m1: estimate beta
    if (use_weight == 1) {
      beta = wpanel_beta(XX, xxinv, W, YY, FE); //  xxinv is xwxinv
    } else {
      beta = panel_beta(XX, xxinv, YY, FE);
    }

    covar_fit.zeros();
    for (int i = 0; i < p; i++) {
      covar_fit = covar_fit + XX.slice(i) * beta(i);
    }

    // m2: estimate interactive fe, additive fe, and mu
    // To simplify, we ignore extra FE in this step for now
    if (use_weight == 1) {
      U = wE_adj(YY - covar_fit, FE, W, I); // e step: expeactation
    } else {
      U = E_adj(YY - covar_fit, FE, I);
    }

    if (use_weight == 1 && stop_burnin == 0) {
      r_burnin = d - niter;
      if (r_burnin <= r) {
        r_burnin = r;
      }
      ife_inner = ife(U, force, 0, r_burnin, 0, 0);
    } else {
      ife_inner = ife(U, force, 0, r, 0, 0);
    }

    FE = as<arma::mat>(ife_inner["FE"]);
    fit = covar_fit + FE; // overall fe */
    // fit_FE_inter_use = as<arma::mat>(ife_inner["FE_inter_use"]);
    if (use_weight == 1) {
      dif = arma::norm(W % (fit - fit_old), "fro") /
            arma::norm(W % (fit_old), "fro");
    } else {
      dif = arma::norm(fit - fit_old, "fro") / arma::norm(fit_old, "fro");
    }
    // dif_FE_inter_use = arma::norm(fit_FE_inter_use - fit_FE_inter_use_old, "fro") / arma::norm(fit_FE_inter_use_old, "fro");
    // fit_FE_inter_use_old = fit_FE_inter_use;
    fit_old = fit;
    niter = niter + 1;

    if (dif <= tolerate && niter <= d && use_weight == 1 && stop_burnin == 0) {
      stop_burnin = 1;
      dif = 1.0;
      niter = 0;
      fit = Y0;
      fit_old = fit;
    }
  }
  e = FE_adj(Y - fit, I);

  FE_inter_use = as<arma::mat>(ife_inner["FE_inter_use"]);
  if (arma::accu(abs(FE_inter_use)) < 1e-10) {
    validF = 0;
  }
  List result;
  mu = as<double>(ife_inner["mu"]);
  result["mu"] = mu;
  result["niter"] = niter;
  result["e"] = e;
  result["beta"] = beta;
  result["fit"] = fit;
  result["validF"] = validF;

  if (force == 1 || force == 3) {
    alpha = as<arma::mat>(ife_inner["alpha"]);
    result["alpha"] = alpha;
  }
  if (force == 2 || force == 3) {
    xi = as<arma::mat>(ife_inner["xi"]);
    result["xi"] = xi;
  }
  L = as<arma::mat>(ife_inner["lambda"]);
  F = as<arma::mat>(ife_inner["factor"]);
  VNT = as<arma::mat>(ife_inner["VNT"]);
  result["lambda"] = L;
  result["factor"] = F;
  result["VNT"] = VNT;
  
  if (use_weight == 1) {
    result["burn_in"] = abs(1 - stop_burnin);
  }
  result["gamma"] = gamma;
  result["kappa"] = kappa;
  return (result);
}
