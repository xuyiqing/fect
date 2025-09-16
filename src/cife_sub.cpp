#include "fect.h"

// [[Rcpp::export]]
arma::mat YY_adj(arma::mat YYYY, arma::mat EEE, arma::mat I, int use_weight,
                 arma::mat W) {
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
List Demean(arma::mat E, int force, arma::cube X_extra_FE,
            std::vector<std::vector<arma::uvec>> extra_FE_index_cache) {
  int T = E.n_rows;
  int N = E.n_cols;
  int p_extra_FE = X_extra_FE.n_slices;
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

  if (p_extra_FE > 0) {
    for (int i = 0; i < p_extra_FE; i++) {
      arma::mat extra_FE_mean(T, N, arma::fill::zeros);
      for (unsigned int g = 0; g < extra_FE_index_cache[i].size(); g++) {
        arma::uvec idx = extra_FE_index_cache[i][g];
        if (idx.n_elem == 0)
          continue;
        double m = arma::mean(E.elem(idx));
        extra_FE_mean.elem(idx).fill(m);
      }
      fit = fit + extra_FE_mean;
      E = E - extra_FE_mean;
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
List fixed_effects_part(
    arma::mat E, int force, arma::cube X_extra_FE,
    std::vector<std::vector<arma::uvec>> extra_FE_index_cache) {
  List result;
  result = Demean(E, force, X_extra_FE, extra_FE_index_cache);
  return (result);
}

// [[Rcpp::export]]
arma::mat Gamma(const arma::mat E, const arma::mat Z, const arma::mat zzinv,
                arma::uvec gamma_t_group) {
  const int T = E.n_rows;
  const int N = E.n_cols;
  const int p = Z.n_cols;
  arma::uvec uniq = arma::unique(gamma_t_group);
  const int G = uniq.n_elem;

  arma::mat sumE(G, N, arma::fill::zeros);
  arma::ivec cnt(G, arma::fill::zeros);
  for (int t = 0; t < T; ++t) {
    int g = static_cast<int>(gamma_t_group(t));
    sumE.row(g) += E.row(t);
    cnt(g)++;
  }

  arma::mat gamma_group(p, G, arma::fill::zeros); // p x G
  for (int g = 0; g < G; ++g) {
    if (cnt(g) == 0)
      continue;
    arma::vec y_bar = (sumE.row(g) / double(cnt(g))).t(); // N x 1
    gamma_group.col(g) = zzinv * (Z.t() * y_bar);         // p x 1
  }

  arma::mat gamma(p, T, arma::fill::zeros);
  for (int t = 0; t < T; ++t) {
    int g = static_cast<int>(gamma_t_group(t));
    gamma.col(t) = gamma_group.col(g);
  }
  return gamma; // p x T
}

// [[Rcpp::export]]
List gamma_part(arma::mat E, arma::mat Z, arma::mat zzinv,
                arma::uvec gamma_t_group) {
  arma::mat gamma = Gamma(E, Z, zzinv, gamma_t_group);

  arma::mat fit = (Z * gamma).t();

  E = E - fit;

  List result;
  result["E"] = E;
  result["fit"] = fit;
  result["gamma"] = gamma.t();
  return (result);
}

// [[Rcpp::export]]
arma::mat Kappa(const arma::mat E, const arma::mat Q, const arma::mat qqinv,
                arma::uvec kappa_i_group) {
  const int T = E.n_rows;
  const int N = E.n_cols;
  const int q = Q.n_rows;
  const arma::uvec uniq = arma::unique(kappa_i_group);
  const int G = uniq.n_elem;

  arma::mat sumE(T, G, arma::fill::zeros);
  arma::ivec cnt(G, arma::fill::zeros);
  for (int i = 0; i < N; ++i) {
    int g = static_cast<int>(kappa_i_group(i));
    sumE.col(g) += E.col(i);
    cnt(g)++;
  }

  arma::mat kappa_group(G, q, arma::fill::zeros);
  for (int g = 0; g < G; ++g) {
    if (cnt(g) == 0)
      continue;
    arma::vec ybar = sumE.col(g) / double(cnt(g));   // T x 1
    kappa_group.row(g) = (ybar.t() * Q.t()) * qqinv; // 1 x q
  }

  arma::mat kappa(N, q, arma::fill::zeros);
  for (int i = 0; i < N; ++i) {
    int g = static_cast<int>(kappa_i_group(i));
    kappa.row(i) = kappa_group.row(g);
  }
  return kappa;
}

// [[Rcpp::export]]
List kappa_part(arma::mat E, arma::mat Q, arma::mat qqinv,
                arma::uvec kappa_i_group) {
  arma::mat kappa = Kappa(E, Q, qqinv, kappa_i_group);

  arma::mat fit = (kappa * Q).t();

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
List cife_iter(arma::cube XX, arma::mat xxinv, arma::cube X_extra_FE,
               arma::cube X_Z, arma::cube X_Q, arma::cube X_gamma,
               arma::cube X_kappa, Rcpp::List Zgamma_id, Rcpp::List kappaQ_id,
               arma::mat Y, arma::mat Y0, arma::mat I, arma::mat W,
               arma::mat beta0, int force, int r, double tolerate,
               int max_iter = 1000) {
  int T = Y.n_rows;
  int N = Y.n_cols;
  int p = XX.n_slices;
  int p_gamma = X_gamma.n_slices;
  int p_kappa = X_kappa.n_slices;
  double dif = 1.0;
  double dif1 = 1.0;
  double dif2 = 1.0;
  double dif2_tol = 0.0;
  double dif3 = 1.0;
  double dif3_tol = 0.0;
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
  std::vector<List> result2(p_gamma);
  std::vector<List> result3(p_kappa);
  List result4;
  List result5;

  arma::mat fit(T, N, arma::fill::zeros);
  arma::mat fit_old(T, N, arma::fill::ones);

  arma::mat fit1(T, N, arma::fill::zeros);
  arma::mat fit1_old(T, N, arma::fill::ones);

  arma::cube fit2(T, N, p_gamma, arma::fill::zeros);
  arma::mat fit2sum = arma::sum(fit2, 2);
  arma::mat fit2leftsum(T, N, arma::fill::zeros);
  arma::cube fit2_old(T, N, p_gamma, arma::fill::ones);

  arma::cube fit3(T, N, p_kappa, arma::fill::zeros);
  arma::mat fit3sum = arma::sum(fit3, 2);
  arma::mat fit3leftsum(T, N, arma::fill::zeros);
  arma::cube fit3_old(T, N, p_kappa, arma::fill::ones);

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

  FE = fit - fit1 - fit2sum - fit3sum; // fit4 + fit5

  int N_time_inv = X_Z.n_cols;
  int p_time_inv = X_Z.n_slices;
  arma::mat Z(N_time_inv, p_time_inv, arma::fill::zeros);
  for (int k = 0; k < p_time_inv; ++k) {
    Z.col(k) = X_Z.slice(k).row(0).t();
  }

  std::vector<arma::uvec> gamma_t_group(p_gamma);
  std::vector<arma::uvec> Zgamma_id_raw(p_gamma);
  std::vector<arma::mat> zzinv(p_gamma);
  for (int k = 0; k < p_gamma; ++k) {
    arma::uvec raw_labels(T);
    for (int t = 0; t < T; ++t) {
      raw_labels(t) = static_cast<unsigned int>(X_gamma(t, 0, k));
    }

    arma::uvec uniq = arma::unique(raw_labels);
    arma::uvec mapped(T);
    for (int t = 0; t < T; ++t) {
      arma::uvec pos = arma::find(uniq == raw_labels(t), 1);
      unsigned int idx =
          (pos.n_elem > 0) ? static_cast<unsigned int>(pos(0)) : 0u;
      mapped(t) = idx;
    }
    gamma_t_group[k] = std::move(mapped);

    Rcpp::IntegerVector v = Zgamma_id[k];
    arma::uvec raw = Rcpp::as<arma::uvec>(v);
    raw.transform([](arma::uword x) { return x > 0 ? x - 1 : 0; });
    Zgamma_id_raw[k] = std::move(raw);

    zzinv[k] = arma::inv_sympd(Z.cols(Zgamma_id_raw[k]).t() *
                               Z.cols(Zgamma_id_raw[k]));
  }

  int T_time_trend = X_Q.n_rows;
  int p_time_trend = X_Q.n_slices;
  arma::mat Q(p_time_trend, T_time_trend, arma::fill::zeros);
  for (int k = 0; k < p_time_trend; ++k) {
    Q.row(k) = X_Q.slice(k).col(0).t();
  }

  std::vector<arma::uvec> kappa_i_group(p_kappa);
  std::vector<arma::uvec> kappaQ_id_raw(p_kappa);
  std::vector<arma::mat> qqinv(p_kappa);
  for (int k = 0; k < p_kappa; ++k) {
    arma::uvec raw_labels(N);
    for (int i = 0; i < N; ++i) {
      raw_labels(i) = static_cast<unsigned int>(X_kappa(0, i, k));
    }

    arma::uvec uniq = arma::unique(raw_labels);
    arma::uvec mapped(N);
    for (int i = 0; i < N; ++i) {
      arma::uvec pos = arma::find(uniq == raw_labels(i), 1);
      unsigned int idx =
          (pos.n_elem > 0) ? static_cast<unsigned int>(pos(0)) : 0u;
      mapped(i) = idx;
    }
    kappa_i_group[k] = std::move(mapped);

    Rcpp::IntegerVector v = kappaQ_id[k];
    arma::uvec raw = Rcpp::as<arma::uvec>(v);
    raw.transform([](arma::uword x) { return x > 0 ? x - 1 : 0; });
    kappaQ_id_raw[k] = std::move(raw);

    // qqinv should be (Q Q^T)^{-1} with Q having shape (q x T)
    qqinv[k] = arma::inv_sympd(Q.rows(kappaQ_id_raw[k]) *
                               Q.rows(kappaQ_id_raw[k]).t());
  }

  int p_extra_FE = X_extra_FE.n_slices;
  std::vector<std::vector<arma::uvec>> extra_FE_index_cache;
  extra_FE_index_cache.resize(p_extra_FE);

  if (p_extra_FE > 0) {
    for (int i = 0; i < p_extra_FE; i++) {
      arma::mat lab = X_extra_FE.slice(i);
      arma::vec uniq = arma::unique(arma::vectorise(lab));
      std::vector<arma::uvec> idx_list;
      idx_list.reserve(uniq.n_elem);
      for (unsigned int g = 0; g < uniq.n_elem; g++) {
        arma::uvec idx = arma::find(lab == uniq(g));
        if (idx.n_elem > 0) {
          idx_list.push_back(idx);
        }
      }
      extra_FE_index_cache[i] = idx_list;
    }
  }

  int stop_burnin = 0;
  while (((dif > tolerate) || (dif1 > tolerate) || dif2_tol || dif3_tol ||
          (dif4 > tolerate) || (dif5 > tolerate)) &&
         (niter <= max_iter)) {
    YYY = E_adj(Y, fit, I);
    result1 = beta_part(YYY - fit2sum - fit3sum - FE, XX, xxinv, W, use_weight);
    fit1 = as<arma::mat>(result1["fit"]);

    for (int k = 0; k < p_gamma; ++k) {
      fit2leftsum = fit2sum - fit2.slice(k);
      YY = YY_adj(YYY - fit1 - fit3sum - FE - fit2leftsum, fit2.slice(k), I,
                  use_weight, W);
      result2[k] =
          gamma_part(YY, Z.cols(Zgamma_id_raw[k]), zzinv[k], gamma_t_group[k]);
      fit2.slice(k) = as<arma::mat>(result2[k]["fit"]);
      fit2sum = arma::sum(fit2, 2);
    }

    for (int k = 0; k < p_kappa; ++k) {
      fit3leftsum = fit3sum - fit3.slice(k);
      YY = YY_adj(YYY - fit1 - fit3leftsum - FE - fit2sum, fit3.slice(k), I,
                  use_weight, W);
      result3[k] =
          kappa_part(YY, Q.rows(kappaQ_id_raw[k]), qqinv[k], kappa_i_group[k]);
      fit3.slice(k) = as<arma::mat>(result3[k]["fit"]);
      fit3sum = arma::sum(fit3, 2);
    }

    YY = YY_adj(YYY - fit1 - fit2sum - fit3sum - fit5, FE - fit5, I, use_weight,
                W);
    result4 = fixed_effects_part(YY, force, X_extra_FE, extra_FE_index_cache);
    fit4 = as<arma::mat>(result4["fit"]);

    YY = YY_adj(YYY - fit1 - fit2sum - fit3sum - fit4, FE - fit4, I, use_weight,
                W);
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
    fit = fit1 + fit2sum + fit3sum + fit4 + fit5;

    if (use_weight == 1) {
      dif = arma::norm(W % (fit - fit_old), "fro") /
            arma::norm(W % (fit_old), "fro");
    } else {
      dif = arma::norm(fit - fit_old, "fro") / arma::norm(fit_old, "fro");
    }

    dif1 = arma::norm(fit1 - fit1_old, "fro") / arma::norm(fit1_old, "fro");

    dif2_tol = 0.0;
    for (int k = 0; k < p_gamma; ++k) {
      dif2 = arma::norm(fit2.slice(k) - fit2_old.slice(k), "fro") /
             arma::norm(fit2_old.slice(k), "fro");
      dif2_tol = dif2_tol || (dif2 > tolerate);
    }

    for (int k = 0; k < p_kappa; ++k) {
      dif3 = arma::norm(fit3.slice(k) - fit3_old.slice(k), "fro") /
             arma::norm(fit3_old.slice(k), "fro");
      dif3_tol = dif3_tol || (dif3 > tolerate);
    }

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

  List gamma_list(p_gamma);
  for (int k = 0; k < p_gamma; ++k)
    gamma_list[k] = result2[k]["gamma"];
  result["gamma"] = gamma_list;

  List kappa_list(p_kappa);
  for (int k = 0; k < p_kappa; ++k)
    kappa_list[k] = result3[k]["kappa"];
  result["kappa"] = kappa_list;
  return (result);
}
