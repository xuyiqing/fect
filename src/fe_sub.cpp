#include "fect.h"

// functions for fixed effects

/* unbalanced panel: response demean function */
// [[Rcpp::export]]
List Y_demean(const arma::mat& Y, int force) {
  int T = Y.n_rows;
  int N = Y.n_cols;
  double mu_Y = 0;
  arma::mat alpha_Y(N, 1, arma::fill::zeros);
  arma::mat xi_Y(T, 1, arma::fill::zeros);
  arma::mat YY = Y;
  mu_Y = accu(YY) / (N * T);
  if (force == 0) {
    YY = YY - mu_Y;
  }
  /* unit fixed effects */
  if (force == 1) {
    alpha_Y = mean(YY, 0).t();           // colMeans, (N * 1) matrix
    YY = YY - repmat(alpha_Y.t(), T, 1); // (T * N) matrix
  }
  /* time fixed effects  */
  if (force == 2) {
    xi_Y = mean(YY, 1); // rowMeans, (N * 1) matrix
    YY = YY - repmat(xi_Y, 1, N);
  }
  if (force == 3) {
    alpha_Y = mean(YY, 0).t();
    xi_Y = mean(YY, 1);
    YY = YY - repmat(alpha_Y.t(), T, 1) - repmat(xi_Y, 1, N) + mu_Y;
  }

  List result;

  result["mu_Y"] = mu_Y;
  result["YY"] = YY;

  if (force == 1 || force == 3) {
    result["alpha_Y"] = alpha_Y;
  }
  if (force == 2 || force == 3) {
    result["xi_Y"] = xi_Y;
  }
  return (result);
}

/* Weighted Demean Outcome Matrix*/
/* return fixed effects */
// [[Rcpp::export]]
List Y_wdemean(const arma::mat& Y, const arma::mat& W, int force) {
  int T = Y.n_rows;
  int N = Y.n_cols;
  double mu_Y = 0;
  arma::mat alpha_Y(N, 1, arma::fill::zeros);
  arma::mat xi_Y(T, 1, arma::fill::zeros);
  arma::mat YY = Y;

  mu_Y = accu(YY % W) / accu(W);
  if (force == 0) {
    YY = YY - mu_Y;
  }

  /* unit fixed effects */
  if (force == 1) {
    alpha_Y = (sum(YY % W, 0) / sum(W, 0)).t(); // colMeans, (N * 1) matrix
    YY = YY - repmat(alpha_Y.t(), T, 1);        // (T * N) matrix
  }

  /* time fixed effects  */
  if (force == 2) {
    xi_Y = (sum(YY % W, 1) / sum(W, 1)); // rowMeans, (N * 1) matrix
    YY = YY - repmat(xi_Y, 1, N);
  }

  if (force == 3) {
    alpha_Y = (sum(YY % W, 0) / sum(W, 0)).t();
    xi_Y = (sum(YY % W, 1) / sum(W, 1));
    YY = YY - repmat(alpha_Y.t(), T, 1) - repmat(xi_Y, 1, N) + mu_Y;
  }

  List result;
  result["mu_Y"] = mu_Y;
  result["YY"] = YY;
  if (force == 1 || force == 3) {
    result["alpha_Y"] = alpha_Y;
  }
  if (force == 2 || force == 3) {
    result["xi_Y"] = xi_Y;
  }
  return (result);
}

/* estimate additive fe for unbalanced panel, without covariates */
// [[Rcpp::export]]
List fe_add(const arma::mat& alpha_Y, const arma::mat& xi_Y, double mu_Y, int T, int N,
            int force) {
  arma::mat FE_ad(T, N, arma::fill::zeros);
  double mu = 0;
  arma::mat alpha(N, 1, arma::fill::zeros);
  arma::mat xi(T, 1, arma::fill::zeros);

  mu = mu_Y;
  if (force == 1 || force == 3) {
    alpha = alpha_Y - mu_Y;
  }
  if (force == 2 || force == 3) {
    xi = xi_Y - mu_Y;
  }

  FE_ad = FE_ad + mu;

  if (force == 1 || force == 3) {
    FE_ad = FE_ad + repmat(alpha.t(), T, 1);
  }
  if (force == 2 || force == 3) {
    FE_ad = FE_ad + repmat(xi, 1, N);
  }

  List result;
  result["mu"] = mu;
  result["FE_ad"] = FE_ad;
  if (force == 1 || force == 3) {
    result["alpha"] = alpha;
  }
  if (force == 2 || force == 3) {
    result["xi"] = xi;
  }

  return (result);
}

/* Obtain factors and loading given error */
// [[Rcpp::export]]
List panel_factor(const arma::mat& E, int r) {
  int T = E.n_rows;
  int N = E.n_cols;
  arma::mat factor(T, r, arma::fill::zeros);
  arma::mat lambda(N, r, arma::fill::zeros);
  arma::mat FE(T, N, arma::fill::zeros);
  arma::mat VNT(r, r, arma::fill::zeros);
  arma::mat U;
  arma::vec s;
  arma::mat V;
  if (T < N) {
    /*arma::mat EE = E * E.t() /(N * T) ;
    arma::eig_sym(s, U, EE) ;
    factor = U.tail_cols(r) * sqrt(double(T)) ;
    lambda = E.t() * factor/T ;
    VNT = diagmat(s.tail_rows(r)) ;*/
    arma::mat EE = E * E.t() / (N * T);
    arma::svd(U, s, V, EE);
    factor = U.head_cols(r) * sqrt(double(T));
    lambda = E.t() * factor / T;
    VNT = diagmat(s.head_rows(r));
  } else {
    /*arma::mat EE = E.t() * E / (N * T) ;
    arma::eig_sym(s, U, EE) ;
    lambda = U.tail_cols(r) * sqrt(double(N)) ;
    factor = E * lambda / N ;
    VNT = diagmat(s.tail_rows(r)) ;*/
    arma::mat EE = E.t() * E / (N * T);
    svd(U, s, V, EE);
    lambda = U.head_cols(r) * sqrt(double(N));
    factor = E * lambda / N;
    VNT = diagmat(s.head_rows(r));
  }
  FE = factor * lambda.t();
  List result;
  result["lambda"] = lambda;
  result["factor"] = factor;
  result["VNT"] = VNT;
  result["FE"] = FE;
  return (result);
}

/* Obtain interactive fe directly */
// [[Rcpp::export]]
arma::mat panel_FE(const arma::mat& E, double lambda, int hard) {
  int T = E.n_rows;
  int N = E.n_cols;
  int r = T;
  if (T >= N) {
    r = N;
  }

  arma::mat FE(T, N, arma::fill::zeros);
  arma::mat D(r, r, arma::fill::zeros);
  arma::mat U;
  arma::vec s;
  arma::mat V;
  arma::svd(U, s, V, E / (T * N));

  for (int i = 0; i < r; i++) {
    if (s(i) > lambda) {
      if (hard == 1) {
        D(i, i) = s(i); // hard impute
      } else {
        D(i, i) = s(i) - lambda; // soft impute
      }
    } else {
      D(i, i) = 0;
    }
  }
  if (T >= N) {
    arma::mat UU = U.cols(0, r - 1);
    FE = UU * D * V.t() * (T * N);
  } else {
    arma::mat VV = V.cols(0, r - 1);
    FE = U * D * VV.t() * (T * N);
  }
  return (FE);
}

/* factor analysis: mu add ife*/
// [[Rcpp::export]]
List ife(const arma::mat& E, int force,
         int mc, // whether pac or mc method
         int r, int hard, double lambda) {
  int T = E.n_rows;
  int N = E.n_cols;

  arma::mat VNT(r, r);
  arma::mat EE(T, N, arma::fill::zeros);
  arma::mat FE_add_use(T, N, arma::fill::zeros);
  arma::mat FE_inter_use(T, N, arma::fill::zeros);
  arma::mat FE(T, N, arma::fill::zeros);

  arma::mat F(T, r, arma::fill::zeros);
  arma::mat L(N, r, arma::fill::zeros);

  double mu_E = 0;
  arma::mat alpha_E(N, 1, arma::fill::zeros);
  arma::mat xi_E(T, 1, arma::fill::zeros);

  double mu = 0;
  arma::mat alpha(N, 1, arma::fill::zeros);
  arma::mat xi(T, 1, arma::fill::zeros);

  List E_ad;
  List E_fe_ad;
  List pf;
  E_ad = Y_demean(E, force);
  EE = as<arma::mat>(E_ad["YY"]);
  mu_E = as<double>(E_ad["mu_Y"]);
  if (force == 1 || force == 3) {
    alpha_E = as<arma::mat>(E_ad["alpha_Y"]);
  }
  if (force == 2 || force == 3) {
    xi_E = as<arma::mat>(E_ad["xi_Y"]);
  }
  E_fe_ad = fe_add(alpha_E, xi_E, mu_E, T, N, force);
  FE_add_use = as<arma::mat>(E_fe_ad["FE_ad"]); // additive fe

  if (r > 0) {
    if (mc == 0) {
      pf = panel_factor(EE, r);
      F = as<arma::mat>(pf["factor"]);
      L = as<arma::mat>(pf["lambda"]);
      VNT = as<arma::mat>(pf["VNT"]);
      FE_inter_use = F * L.t(); // interactive fe
    } else {
      FE_inter_use = panel_FE(EE, lambda, hard);
    }
  }

  FE = FE_add_use + FE_inter_use;

  List result;

  mu = as<double>(E_fe_ad["mu"]);
  result["mu"] = mu;
  result["FE"] = FE;
  if (force == 1 || force == 3) {
    alpha = as<arma::mat>(E_fe_ad["alpha"]);
    result["alpha"] = alpha;
  }
  if (force == 2 || force == 3) {
    xi = as<arma::mat>(E_fe_ad["xi"]);
    result["xi"] = xi;
  }
  if (r > 0) {
    if (mc == 0) {
      result["lambda"] = L;
      result["factor"] = F;
      result["VNT"] = VNT;
    }
    result["FE_inter_use"] = FE_inter_use;
  }
  return (result);
}
