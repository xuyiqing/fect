#include "fect.h"

// [[Rcpp::export]]
arma::mat YY_adj(arma::mat YYYY, arma::mat EEE, arma::mat I, int use_weight, arma::mat W) {
  arma::mat YY = YYYY;
  if (use_weight == 1) {
    YY = wE_adj(YYYY, EEE, W, I); // e step: expeactation
  } else {
    YY = E_adj(YYYY, EEE, I);
  }
  return YY;
}

/* unbalanced panel: response demean function */
// [[Rcpp::export]]
List Demean(arma::mat E, int force, arma::cube X_sfe, std::vector<std::vector<arma::uvec>> sfe_index_cache) {
  int T = E.n_rows;
  int N = E.n_cols;
  int p_sfe = X_sfe.n_slices;
  double mu = 0;
  arma::mat alpha_E(N, 1, arma::fill::zeros);
  arma::mat xi_E(T, 1, arma::fill::zeros);
  arma::mat fit(T, N, arma::fill::zeros);

  mu = accu(E) / (N * T);
  if (force == 0) {
    E = E - mu;
  }
  fit = fit + mu;
  /* unit fixed effects */
  if (force == 1) {
    alpha_E = mean(E, 0).t();
    fit = fit + repmat((alpha_E - mu).t(), T, 1);
    E = E - repmat(alpha_E.t(), T, 1);
  }
  /* time fixed effects  */
  if (force == 2) {
    xi_E = mean(E, 1);
    fit = fit + repmat((xi_E - mu), 1, N);
    E = E - repmat(xi_E, 1, N);
  }
  if (force == 3) {
    alpha_E = mean(E, 0).t();
    fit = fit + repmat((alpha_E - mu).t(), T, 1);

    xi_E = mean(E, 1);
    fit = fit + repmat((xi_E - mu), 1, N);
    E = E - repmat(alpha_E.t(), T, 1) - repmat(xi_E, 1, N) + mu;
  }

  if (p_sfe > 0) {
    for (int i = 0; i < p_sfe; i++) {
      arma::mat sfe_mean(T, N, arma::fill::zeros);
      for (unsigned int g = 0; g < sfe_index_cache[i].size(); g++) {
        arma::uvec idx = sfe_index_cache[i][g];
        if (idx.n_elem == 0)
          continue;
        double m = arma::mean(E.elem(idx));
        sfe_mean.elem(idx).fill(m);
      }
      fit = fit + sfe_mean;
      E = E - sfe_mean;
    }
  }

  List result;

  if (force == 1 || force == 3) {
    result["alpha"] = alpha_E - mu;
  }
  if (force == 2 || force == 3) {
    result["xi"] = xi_E - mu;
  }

  result["E"] = E;
  result["mu"] = mu;
  result["fit"] = fit;
  return (result);
}

// [[Rcpp::export]]
List fixed_effects_part(arma::mat E, int force, arma::cube X_sfe, std::vector<std::vector<arma::uvec>> sfe_index_cache) {
  List result;
  result = Demean(E, force, X_sfe, sfe_index_cache);
  return (result);
}

// [[Rcpp::export]]
arma::mat Gamma(const arma::mat E,const arma::mat Z, const arma::mat zzinv) {
  const int T = E.n_rows;
  arma::mat gamma(Z.n_cols, T, arma::fill::zeros); // p x T
  for (int t = 0; t < T; ++t) {
    arma::vec y_t = E.row(t).t();                 // N x 1
    gamma.col(t)  = zzinv * (Z.t() * y_t);       // p x 1
  }
  return gamma; // p x T  return zeta;
}

// [[Rcpp::export]]
List gamma_part(arma::mat E, arma::mat Z, arma::mat zzinv) {
  arma::mat gamma = Gamma(E, Z, zzinv);

  arma::mat fit = (Z * gamma).t();

  E = E - fit;

  List result;
  result["E"] = E;
  result["fit"] = fit;
  result["gamma"] = gamma.t();
  return (result);
}

// [[Rcpp::export]]
arma::mat Kappa(const arma::mat E,const arma::mat F, const arma::mat ffinv) {
  arma::mat kappa = E.t() * F.t() * ffinv;
  return kappa; // N x p_time_trend
}

// [[Rcpp::export]]
List kappa_part(arma::mat E, arma::mat F, arma::mat ffinv) {
  arma::mat kappa = Kappa(E, F, ffinv);

  arma::mat fit = (kappa * F).t();

  E = E - fit;

  List result;
  result["E"] = E;
  result["fit"] = fit;
  result["kappa"] = kappa.t();
  return (result);
}

/* Obtain beta */
// [[Rcpp::export]]
arma::mat Beta(arma::cube X, arma::mat xxinv, arma::mat E) {
  int p = X.n_slices;
  arma::mat xy(p, 1, arma::fill::zeros);
  for (int k = 0; k < p; k++) {
    xy(k) = trace(crossprod(X.slice(k), E));
  }
  return (xxinv * xy);
}

/* Obtain beta given weights */
// [[Rcpp::export]]
arma::mat WBeta(arma::cube X, arma::mat xwxinv, arma::mat w, arma::mat E) {
  int p = X.n_slices;
  arma::mat xwy(p, 1, arma::fill::zeros);
  arma::mat w_sr = pow(w, 0.5);
  for (int k = 0; k < p; k++) {
    xwy(k) = trace(crossprod(w_sr % X.slice(k), w_sr % E));
  }
  return (xwxinv * xwy);
}
// [[Rcpp::export]]
List beta_part(arma::mat E, arma::cube XX, arma::mat xxinv, arma::mat W,
               int use_weight) {
  int T = E.n_rows;
  int N = E.n_cols;
  int p = XX.n_slices;
  arma::mat beta(p, 1, arma::fill::zeros);
  if (use_weight == 1) {
    beta = WBeta(XX, xxinv, W, E);
  } else {
    beta = Beta(XX, xxinv, E);
  }

  arma::mat fit(T, N, arma::fill::zeros);
  for (int i = 0; i < p; i++) {
    fit = fit + XX.slice(i) * beta(i);
  }

  E = E - fit;

  List result;
  result["E"] = E;
  result["fit"] = fit;
  result["beta"] = beta;
  return (result);
}

// [[Rcpp::export]]
List ife_part(arma::mat E, int r) {
  int T = E.n_rows;
  int N = E.n_cols;

  List pf;
  arma::mat F(T, r, arma::fill::zeros);
  arma::mat L(N, r, arma::fill::zeros);
  arma::mat VNT(r, r);
  pf = panel_factor(E, r);
  F = as<arma::mat>(pf["factor"]);
  L = as<arma::mat>(pf["lambda"]);
  VNT = as<arma::mat>(pf["VNT"]);

  arma::mat fit(T, N, arma::fill::zeros);
  fit = F * L.t();
  E = E - fit;
  List result;
  result["factor"] = F;
  result["lambda"] = L;
  result["VNT"] = VNT;
  result["fit"] = fit;
  result["E"] = E;
  return (result);
}

/* Obtain cife; */
// [[Rcpp::export]]
List cife_iter(arma::cube XX, arma::mat xxinv, arma::cube X_sfe,
               arma::cube X_time_inv,
               arma::cube X_time_trend, 
               arma::mat Y,
               arma::mat Y0, arma::mat I, arma::mat W, arma::mat beta0,
               int force, int r, double tolerate, int max_iter = 1000) {
  int T = Y.n_rows;
  int N = Y.n_cols;
  int p = XX.n_slices;
  double dif = 1.0;
  double dif1 = 1.0;
  double dif2 = 1.0;
  double dif3 = 1.0;
  double dif4 = 1.0;
  double dif5 = 1.0;
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

  arma::mat e(T, N, arma::fill::zeros); // residual
  arma::mat YY = Y;
  arma::mat YYY = Y;

  List result1;
  List result2;
  List result3;
  List result4;
  List result5;
  arma::mat fit(T, N, arma::fill::zeros);
  arma::mat fit_old(T, N, arma::fill::ones);
  arma::mat fit1(T, N, arma::fill::zeros);
  arma::mat fit1_old(T, N, arma::fill::ones);
  arma::mat fit2(T, N, arma::fill::zeros);
  arma::mat fit2_old(T, N, arma::fill::ones);
  arma::mat fit3(T, N, arma::fill::zeros);
  arma::mat fit3_old(T, N, arma::fill::ones);
  arma::mat fit4(T, N, arma::fill::zeros);
  arma::mat fit4_old(T, N, arma::fill::ones);
  arma::mat fit5(T, N, arma::fill::zeros);
  arma::mat fit5_old(T, N, arma::fill::ones);
  arma::mat FE(T, N, arma::fill::zeros);

  fit = Y0;
  fit_old = fit;
  for (int i = 0; i < p; i++) {
    fit1 = fit1 + XX.slice(i) * beta0(i);
  }

  FE = fit - fit1 - fit2 - fit3; // fit4 + fit5

  int N_time_inv = X_time_inv.n_cols;
  int p_time_inv = X_time_inv.n_slices;
  arma::mat Z(N_time_inv, p_time_inv, arma::fill::zeros);
  for (int k = 0; k < p_time_inv; ++k) {
    Z.col(k) = X_time_inv.slice(k).row(0).t();
  }
  arma::mat zzinv = arma::inv_sympd(Z.t() * Z);

  const int T_time_trend = X_time_trend.n_rows;
  const int p_time_trend = X_time_trend.n_slices;
  arma::mat F(p_time_trend, T_time_trend, arma::fill::zeros);
  for (int k = 0; k < p_time_trend; ++k) {
    F.row(k) = X_time_trend.slice(k).col(0).t();
  }
  arma::mat ffinv = arma::inv_sympd(F * F.t());

  int p_sfe = X_sfe.n_slices;
  std::vector<std::vector<arma::uvec>> sfe_index_cache;
  sfe_index_cache.resize(p_sfe);

  if (p_sfe > 0) {
    for (int i = 0; i < p_sfe; i++) {
      arma::mat lab = X_sfe.slice(i);
      arma::vec uniq = arma::unique(arma::vectorise(lab));
      std::vector<arma::uvec> idx_list;
      idx_list.reserve(uniq.n_elem);
      for (unsigned int g = 0; g < uniq.n_elem; g++) {
        arma::uvec idx = arma::find(lab == uniq(g));
        if (idx.n_elem > 0) {
          idx_list.push_back(idx);
        }
      }
      sfe_index_cache[i] = idx_list;
    }
  }

  int stop_burnin = 0;
  while ((dif > tolerate || dif1 > tolerate || dif2 > tolerate ||
          dif3 > tolerate || dif4 > tolerate || dif5 > tolerate) &&
         niter <= max_iter) {
    YYY = E_adj(Y, fit, I);
    result1 = beta_part(YYY - fit2 - fit3 - FE, XX, xxinv, W, use_weight);
    fit1 = as<arma::mat>(result1["fit"]);

    YY = YY_adj(YYY - fit1 - fit3 - FE, fit2, I, use_weight, W);
    result2 = gamma_part(YY, Z, zzinv);
    fit2 = as<arma::mat>(result2["fit"]);

    YY = YY_adj(YYY - fit1 - fit2 - FE, fit3, I, use_weight, W);
    result3 = kappa_part(YY, F, ffinv);
    fit3 = as<arma::mat>(result3["fit"]);

    YY = YY_adj(YYY - fit1 - fit2 - fit3 - fit5, FE - fit5, I, use_weight, W);
    result4 = fixed_effects_part(YY, force, X_sfe, sfe_index_cache);
    fit4 = as<arma::mat>(result4["fit"]);

    YY = YY_adj(YYY - fit1 - fit2 - fit3 - fit4, FE - fit4, I, use_weight, W);
    if (use_weight == 1 && stop_burnin == 0) {
      r_burnin = d - niter;
      if (r_burnin <= r) {
        r_burnin = r;
      }
      result5 = ife_part(YY, r_burnin);
    } else {
      result5 = ife_part(YY, r);
    }
    fit5 = as<arma::mat>(result5["fit"]);

    FE = fit4 + fit5;
    fit = fit1 + fit2 + fit3 + fit4 + fit5;

    if (use_weight == 1) {
      dif = arma::norm(W % (fit - fit_old), "fro") /
            arma::norm(W % (fit_old), "fro");
    } else {
      dif = arma::norm(fit - fit_old, "fro") / arma::norm(fit_old, "fro");
    }
    dif1 = arma::norm(fit1 - fit1_old, "fro") / arma::norm(fit1_old, "fro");
    dif2 = arma::norm(fit2 - fit2_old, "fro") / arma::norm(fit2_old, "fro");
    dif3 = arma::norm(fit3 - fit3_old, "fro") / arma::norm(fit3_old, "fro");
    dif4 = arma::norm(fit4 - fit4_old, "fro") / arma::norm(fit4_old, "fro");
    dif5 = arma::norm(fit5 - fit5_old, "fro") / arma::norm(fit5_old, "fro");
    fit_old = fit;
    fit1_old = fit1;
    fit2_old = fit2;
    fit3_old = fit3;
    fit4_old = fit4;
    fit5_old = fit5;

    niter = niter + 1;
  }
  e = FE_adj(Y - fit, I);
  if (arma::accu(abs(fit5)) < 1e-10) {
    validF = 0;
  }

  List result;
  result["mu"] = result4["mu"];
  result["alpha"] = result4["alpha"];
  result["xi"] = result4["xi"];
  result["niter"] = niter;
  result["e"] = e;
  result["beta"] = result1["beta"];
  result["fit"] = fit;
  result["validF"] = validF;

  result["lambda"] = result5["lambda"];
  result["factor"] = result5["factor"];
  result["VNT"] = result5["VNT"];
  if (use_weight == 1) {
    result["burn_in"] = abs(1 - stop_burnin);
  }
  result["time_invariant"] = result2["gamma"];
  result["time_trend"] = result3["kappa"];
  return (result);
}
