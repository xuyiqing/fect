#include "fect.h"

// core cfe model function

/* Complex Fixed Effects: ub */
// [[Rcpp::export]]
List complex_fe_ub(
    const arma::mat& Y, const arma::mat& Y0, const arma::cube& X_covariates, const arma::cube& X_extra_FE,
    const arma::cube& X_Z, const arma::cube& X_Q, const arma::cube& X_gamma, const arma::cube& X_kappa,
    Rcpp::List Zgamma_id, Rcpp::List kappaQ_id, const arma::mat& I, const arma::mat& W_in,
    const arma::mat& beta0,
    int r, // r > 0, the outcome has a factor-type fixed effect; r = 0 else
    int force, double tol = 1e-5, int max_iter = 1000) {

  arma::mat W = W_in;
  /* Dimensions */
  int T = Y.n_rows;
  int N = Y.n_cols;
  int p = X_covariates.n_slices;
  double obs = accu(I);
  int niter = 0;
  arma::mat factor;
  arma::mat lambda;
  arma::mat VNT;
  arma::mat beta;
  arma::mat U;
  double mu = 0;
  double mu_Y = 0;
  arma::mat alpha(N, 1, arma::fill::zeros);
  arma::mat xi(T, 1, arma::fill::zeros);
  arma::mat fit(T, N, arma::fill::zeros);
  double sigma2 = 0;
  double IC = 0;
  arma::mat invXX;
  arma::mat xtimeinvinv;
  arma::mat xtimetrendinv;
  int use_weight;
  arma::mat WI;
  int burn_in = 0;

  if (Y.n_rows == W.n_rows && Y.n_cols == W.n_cols) {
    use_weight = 1;
    W = W / W.max();
    WI = W % I;
  } else {
    use_weight = 0;
  }

  /* duplicate data */
  arma::mat YY = Y;
  arma::cube XX = X_covariates;
  /* check if XX has enough variation */
  int p1 = p;
  arma::mat X_invar(p, 1, arma::fill::zeros); // =1 if invar

  int flag_var;
  int j = 0;
  for (int i = 0; i < p1; i++) {
    if (use_weight == 1) {
      flag_var = (arma::accu(abs(XX.slice(i))) < 1e-5);
    } else {
      flag_var = (arma::accu(abs(XX.slice(i))) < 1e-5);
    }
    if (flag_var) {
      XX.shed_slice(i);
      X_invar(j, 0) = 1;
      i--;
      p1--;
    }
    j++;
  }

  int validX = 1;
  if (p1 == 0) {
    validX = 0;
    if (force == 0 && r == 0) { // no covariate and force == 0 and r == 0
      if (use_weight == 1) {
        mu_Y = accu(YY % WI) / accu(WI);
        mu = mu_Y;
        YY = FE_adj(YY - mu_Y, I);
      } else {
        mu_Y = accu(YY) / obs;
        mu = mu_Y;
        YY = FE_adj(YY - mu_Y, I);
      }
    }
  }
  /* Main Algorithm */

  /* starting value:  the OLS estimator */
  if (use_weight == 1) {
    invXX = wXXinv(XX, W);
  } else {
    invXX = XXinv(XX); // compute (X'X)^{-1}, outside beta iteration
  }
  List cfe =
      cfe_iter(XX, invXX, X_extra_FE, X_Z, X_Q, X_gamma, X_kappa, Zgamma_id,
                kappaQ_id, YY, Y0, I, W, beta0, force, r, tol, max_iter);

  mu = as<double>(cfe["mu"]);
  beta = as<arma::mat>(cfe["beta"]);
  U = as<arma::mat>(cfe["e"]);
  fit = as<arma::mat>(cfe["fit"]);
  if (use_weight == 1) {
    burn_in = as<int>(cfe["burn_in"]);
  }
  factor = as<arma::mat>(cfe["factor"]);
  lambda = as<arma::mat>(cfe["lambda"]);
  VNT = as<arma::mat>(cfe["VNT"]);

  if (force == 1 || force == 3) {
    alpha = as<arma::mat>(cfe["alpha"]);
  }
  if (force == 2 || force == 3) {
    xi = as<arma::mat>(cfe["xi"]);
  }
  niter = as<int>(cfe["niter"]);

  /* sigma2 and IC */
  // number of estimated parameters
  double np = r * (N + T) - pow(double(r), 2) + p1 + 1;
  if (force == 1) {
    np = np + (N - 1) - r;
  } else if (force == 2) {
    np = np + (T - 1) - r;
  } else if (force == 3) {
    np = np + (N - 1) + (T - 1) - 2 * r;
  }

  sigma2 = trace(U * U.t()) / (obs - np);

  IC = log(sigma2) + np * log(obs) / (obs);

  // PC criterion in Li 2018
  // mean squared error
  double mse = trace(U * U.t()) / obs;

  double m1 = 0;
  if (N < 60) {
    m1 = 60 - N;
  }
  double m2 = 0;
  if (T < 60) {
    m2 = 60 - T;
  }

  double C = (N + m1) * (T + m2) / obs;
  double PC = mse + r * sigma2 * C * (N + T) / (N * T) *
                        log(double(N * T) / double(N + T));
  //-------------------------------#
  // Storage
  //-------------------------------#

  List output;

  if (p > 0) {
    arma::mat beta_total(p, 1);

    if (p > p1) {
      int j4 = 0;
      for (int i = 0; i < p; i++) {
        if (X_invar(i, 0) == 1) {
          beta_total(i, 0) = arma::datum::nan;
        } else {
          beta_total(i, 0) = beta(j4, 0);
          j4++;
        }
      }
    } else {
      beta_total = beta;
    }
    output["beta"] = beta_total;
  }

  output["mu"] = mu;
  output["fit"] = fit;

  if (!(force == 0 && r == 0 && p1 == 0)) {
    output["niter"] = niter;
  }
  if (force == 1 || force == 3) {
    output["alpha"] = alpha;
  }
  if (force == 2 || force == 3) {
    output["xi"] = xi;
  }
  if (r > 0) {
    output["factor"] = factor;
    output["lambda"] = lambda;
    output["VNT"] = VNT;
    if (use_weight == 1) {
      output["burn_in"] = burn_in;
    }
  }
  output["residuals"] = U;
  output["sigma2"] = sigma2;
  output["IC"] = IC;
  output["PC"] = PC;
  output["validX"] = validX;
  output["gamma"] = cfe["gamma"];
  output["kappa"] = cfe["kappa"];
  return (output);
}
