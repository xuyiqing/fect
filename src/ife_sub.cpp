#include "fect.h"

// sub-functions for unbalanced ife model

/* Obtain additive fe for ub data; assume r=0, without covar */
// [[Rcpp::export]]
List fe_ad_iter(const arma::mat& Y, const arma::mat& Y0, const arma::mat& I, const arma::mat& W, int force,
                double tolerate, int max_iter = 500) {

  int T = Y.n_rows;
  int N = Y.n_cols;
  double mu = 0;
  double dif = 1.0;
  int niter = 0;
  int use_weight;

  if (Y.n_rows == W.n_rows && Y.n_cols == W.n_cols) {
    use_weight = 1;
  } else {
    use_weight = 0;
  }

  arma::mat fit = Y0; // initial value
  arma::mat fit_old = Y0;
  arma::mat alpha(N, 1, arma::fill::zeros);
  arma::mat xi(T, 1, arma::fill::zeros);
  double mu_Y = 0;
  arma::mat alpha_Y(N, 1, arma::fill::zeros);
  arma::mat xi_Y(T, 1, arma::fill::zeros);

  arma::mat e(T, N, arma::fill::zeros); // residual
  arma::mat YY = Y;

  List Y_ad;
  List Y_fe_ad;

  while (dif > tolerate && niter <= 500) {
    if (use_weight == 1) {
      YY = wE_adj(Y, fit, W, I); // e step: expeactation
    } else {
      YY = E_adj(Y, fit, I);
    }

    Y_ad = Y_demean(YY, force);

    mu_Y = as<double>(Y_ad["mu_Y"]);
    if (force == 1 || force == 3) {
      alpha_Y = as<arma::mat>(Y_ad["alpha_Y"]);
    }
    if (force == 2 || force == 3) {
      xi_Y = as<arma::mat>(Y_ad["xi_Y"]);
    }

    Y_fe_ad = fe_add(alpha_Y, xi_Y, mu_Y, T, N, force); // m step: estimate fe

    fit = as<arma::mat>(Y_fe_ad["FE_ad"]);

    if (use_weight == 1) {
      dif = arma::norm(W % (fit - fit_old), "fro") /
            arma::norm(W % (fit_old), "fro");
    } else {
      dif = arma::norm(fit - fit_old, "fro") / arma::norm(fit_old, "fro");
    }

    fit_old = fit;

    niter = niter + 1;
  }
  e = FE_adj(YY - fit, I);

  List result;
  mu = as<double>(Y_fe_ad["mu"]);
  result["mu"] = mu;

  result["fit"] = fit;
  result["niter"] = niter;
  result["e"] = e;
  if (force == 1 || force == 3) {
    alpha = as<arma::mat>(Y_fe_ad["alpha"]);
    result["alpha"] = alpha;
  }
  if (force == 2 || force == 3) {
    xi = as<arma::mat>(Y_fe_ad["xi"]);
    result["xi"] = xi;
  }
  return (result);
}

/* Obtain additive fe for ub data; assume r=0, with covariates */
// [[Rcpp::export]]
List fe_ad_covar_iter(const arma::cube& XX, const arma::mat& xxinv, const arma::mat& Y, const arma::mat& Y0,
                      const arma::mat& I, const arma::mat& beta0, const arma::mat& W, int force,
                      double tolerate, int max_iter = 500) {
  int T = Y.n_rows;
  int N = Y.n_cols;
  int p = XX.n_slices;
  double dif = 1.0;
  int niter = 0;
  int use_weight;

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
  arma::mat fit(T, N, arma::fill::zeros);
  arma::mat fit_old(T, N, arma::fill::ones);

  arma::mat alpha(N, 1, arma::fill::zeros);
  arma::mat xi(T, 1, arma::fill::zeros);
  arma::mat YY = Y;
  arma::mat U(T, N, arma::fill::zeros);
  arma::mat e(T, N, arma::fill::zeros); // residual

  List ife_inner;
  arma::mat covar_fit(T, N, arma::fill::zeros);
  for (int i = 0; i < p; i++) {
    covar_fit = covar_fit + XX.slice(i) * beta(i);
  }

  fit = Y0;
  fit_old = fit;
  arma::mat FE = fit - covar_fit; // initial fixed effects

  while (dif > tolerate && niter <= max_iter) {

    YY = E_adj(Y, fit, I);
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
    if (use_weight == 1) {
      U = wE_adj(YY - covar_fit, FE, W, I); // e step: expeactation
    } else {
      U = E_adj(YY - covar_fit, FE, I);
    }

    ife_inner = ife(U, force, 0, 0, 0, 0);

    FE = as<arma::mat>(ife_inner["FE"]); // overall fe
    fit = covar_fit + FE;

    if (use_weight == 1) {
      dif = arma::norm(W % (fit - fit_old), "fro") /
            arma::norm(W % (fit_old), "fro");
    } else {
      dif = arma::norm(fit - fit_old, "fro") / arma::norm(fit_old, "fro");
    }

    fit_old = fit;
    niter = niter + 1;
  }
  e = FE_adj(YY - fit, I);

  List result;
  mu = as<double>(ife_inner["mu"]);
  result["mu"] = mu;

  result["fit"] = fit;
  result["niter"] = niter;
  result["e"] = e;
  if (p > 0) {
    result["beta"] = beta;
  }
  if (force == 1 || force == 3) {
    alpha = as<arma::mat>(ife_inner["alpha"]);
    result["alpha"] = alpha;
  }
  if (force == 2 || force == 3) {
    xi = as<arma::mat>(ife_inner["xi"]);
    result["xi"] = xi;
  }
  return (result);
}

/* Obtain additive fe for ub data; assume r>0 but p=0*/
// [[Rcpp::export]]
List fe_ad_inter_iter(const arma::mat& Y, const arma::mat& Y0, const arma::mat& I, const arma::mat& W,
                      int force,
                      int mc, // whether pca or mc method
                      int r, int hard, double lambda, double tolerate,
                      int max_iter = 1000) {
  int T = Y.n_rows;
  int N = Y.n_cols;
  double mu = 0;
  double dif = 1.0;
  int niter = 0;
  int validF = 1; // whether has a factor structure
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

  arma::mat VNT(r, r);
  arma::mat FE_inter_use(T, N, arma::fill::zeros);
  arma::mat fit(T, N, arma::fill::zeros);
  arma::mat fit_old(T, N, arma::fill::ones);
  arma::mat e(T, N, arma::fill::zeros); // residual
  arma::mat U(T, N, arma::fill::zeros);

  arma::mat alpha(N, 1, arma::fill::zeros);
  arma::mat xi(T, 1, arma::fill::zeros);

  arma::mat F(T, r, arma::fill::zeros);
  arma::mat L(N, r, arma::fill::zeros);

  arma::mat YY = Y;
  List pf;
  List ife_inner;

  // initial value for ife
  fit = Y0;
  fit_old = fit;
  int stop_burnin = 0;
  while (dif > tolerate && niter <= max_iter) {
    if (use_weight == 1) {
      YY = wE_adj(Y, fit, W, I); // e step: expectation
    } else {
      YY = E_adj(Y, fit, I);
    }
    if (mc == 0) {
      if (use_weight == 1 && stop_burnin == 0) {
        r_burnin = d - niter;
        if (r_burnin <= r) {
          r_burnin = r;
        }
        ife_inner = ife(YY, force, 0, r_burnin, 0, 0);
      } else {
        ife_inner = ife(YY, force, 0, r, 0, 0);
      }
    } else {
      ife_inner = ife(YY, force, 1, 1, hard, lambda);
    }
    fit = as<arma::mat>(ife_inner["FE"]); // new overall fe

    if (use_weight == 1) {
      dif = arma::norm(W % (fit - fit_old), "fro") /
            arma::norm(W % (fit_old), "fro");
    } else {
      dif = arma::norm(fit - fit_old, "fro") / arma::norm(fit_old, "fro");
    }

    fit_old = fit;
    niter = niter + 1;

    if (dif <= tolerate && niter <= d && use_weight == 1 && stop_burnin == 0 &&
        mc == 0) {
      stop_burnin = 1;
      dif = 1.0;
      niter = 0;
      fit = Y0;
      fit_old = fit;
    }
  }
  e = FE_adj(YY - fit, I);
  FE_inter_use = as<arma::mat>(ife_inner["FE_inter_use"]);
  if (arma::accu(abs(FE_inter_use)) < 1e-10) {
    validF = 0;
  }

  List result;
  mu = as<double>(ife_inner["mu"]);
  result["mu"] = mu;
  result["niter"] = niter;
  result["fit"] = fit;
  result["e"] = e;
  result["validF"] = validF;

  if (force == 1 || force == 3) {
    alpha = as<arma::mat>(ife_inner["alpha"]);
    result["alpha"] = alpha;
  }
  if (force == 2 || force == 3) {
    xi = as<arma::mat>(ife_inner["xi"]);
    result["xi"] = xi;
  }
  if (mc == 0) {
    result["lambda"] = as<arma::mat>(ife_inner["lambda"]);
    result["factor"] = as<arma::mat>(ife_inner["factor"]);
    result["VNT"] = as<arma::mat>(ife_inner["VNT"]);
  }
  if (use_weight == 1 && mc == 0) {
    result["burn_in"] = abs(1 - stop_burnin);
  }
  return (result);
}

/* Obtain additive fe for ub data; assume r>0 p>0*/
// [[Rcpp::export]]
List fe_ad_inter_covar_iter(const arma::cube& XX, const arma::mat& xxinv, const arma::mat& Y,
                            const arma::mat& Y0, const arma::mat& I, const arma::mat& W,
                            const arma::mat& beta0, int force,
                            int mc, // whether pca or mc method
                            int r, int hard, double lambda, double tolerate,
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
  // arma::mat fit_FE_inter_use_old(T, N, arma::fill::ones);
  // arma::mat fit_FE_inter_use(T, N, arma::fill::zeros);
  // double dif_FE_inter_use = 1.0;
  // initial value for ife

  // if (hard == 0) {
  //   U = FE_adj(Y - Y0, I) ;
  //   if (mc == 0) {
  //     pf = panel_factor(U, r)  ;
  //     F = as<arma::mat>(pf["factor"]) ;
  //     L = as<arma::mat>(pf["lambda"]) ;
  //     FE_inter_use = F * L.t() ; // interactive fe
  //   }
  //   else {
  //     FE_inter_use = panel_FE(U, lambda, hard) ;
  //   }
  //   fit = Y0 + FE_inter_use ;
  //   fit_old = fit ;
  // } else {
  fit = Y0;
  fit_old = fit;
  //}

  arma::mat covar_fit(T, N, arma::fill::zeros);
  for (int i = 0; i < p; i++) {
    covar_fit = covar_fit + XX.slice(i) * beta(i);
  }
  arma::mat FE = fit - covar_fit;

  int stop_burnin = 0;
  // while ((dif > tolerate || dif_FE_inter_use > tolerate )&& niter <= max_iter) {
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
      ife_inner = ife(U, force, mc, r_burnin, hard, lambda);
    } else {
      ife_inner = ife(U, force, mc, r, hard, lambda);
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

    if (dif <= tolerate && niter <= d && use_weight == 1 && stop_burnin == 0 &&
        mc == 0) {
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
  if (mc == 0) {
    L = as<arma::mat>(ife_inner["lambda"]);
    F = as<arma::mat>(ife_inner["factor"]);
    VNT = as<arma::mat>(ife_inner["VNT"]);
    result["lambda"] = L;
    result["factor"] = F;
    result["VNT"] = VNT;
  }
  if (use_weight == 1 && mc == 0) {
    result["burn_in"] = abs(1 - stop_burnin);
  }
  return (result);
}

/* Main iteration for beta */
// [[Rcpp::export]]
List beta_iter(const arma::cube& X, const arma::mat& xxinv, const arma::mat& Y, int r,
               double tolerate, const arma::mat& beta0, int max_iter) {

  /* beta.new: computed beta under iteration with error precision=tolerate
   factor: estimated factor
   lambda: estimated loadings
   V: the eigenvalues matrix
   e: estimated residuals
   niter: number of interations to achieve convergence */
  int T = Y.n_rows;
  int N = Y.n_cols;
  int p = X.n_slices;
  int b_r = beta0.n_rows;
  double beta_norm = 1.0;
  arma::mat beta(p, 1, arma::fill::zeros);
  if (b_r == p) {
    beta = beta0;
  } // beta should have the same dimension as X, if not it will be reset to 0
  arma::mat beta_old = beta;
  arma::mat VNT(r, r, arma::fill::zeros);
  arma::mat FE(T, N, arma::fill::zeros);

  /* starting value */
  arma::mat U = Y;
  for (int k = 0; k < p; k++) {
    U = U - X.slice(k) * beta(k);
  }
  List pf = panel_factor(U, r);
  arma::mat F = as<arma::mat>(pf["factor"]);
  arma::mat L = as<arma::mat>(pf["lambda"]);

  /* Loop */
  int niter = 0;
  while ((beta_norm > tolerate) && (niter < max_iter)) {
    niter++;
    FE = F * L.t();
    beta = panel_beta(X, xxinv, Y, FE);
    beta_norm = arma::norm(beta - beta_old, "fro");
    beta_old = beta;
    U = Y;
    for (int k = 0; k < p; k++) {
      U = U - X.slice(k) * beta(k);
    }
    pf = panel_factor(U, r);
    F = as<arma::mat>(pf["factor"]);
    L = as<arma::mat>(pf["lambda"]);
  }
  VNT = as<arma::mat>(pf["VNT"]);
  arma::mat e = U - F * L.t();

  /* Storage */
  List result;
  result["niter"] = niter;
  result["beta"] = beta;
  result["e"] = e;
  result["lambda"] = L;
  result["factor"] = F;
  result["VNT"] = VNT;
  return (result);
}
