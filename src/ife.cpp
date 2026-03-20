#include "fect.h"

// core ife model function

/* Interactive Fixed Effects */
// [[Rcpp::export]]
List inter_fe(const arma::mat& Y, const arma::cube& X, int r, int force, const arma::mat& beta0_in,
              double tol = 1e-5, int max_iter = 500) {
  arma::mat beta0 = beta0_in;
  /* Dimensions */
  int b_r = beta0.n_rows;
  int T = Y.n_rows;
  int N = Y.n_cols;
  int p = X.n_slices;
  int obs = T * N;
  int niter = 0;
  arma::mat factor;
  arma::mat lambda;
  arma::mat VNT;
  arma::mat beta;
  arma::mat U;
  double mu = 0;
  double mu_Y = 0;
  arma::mat mu_X(p, 1);
  arma::mat alpha(N, 1, arma::fill::zeros);
  arma::mat alpha_Y(N, 1);
  arma::mat alpha_X(N, p);
  arma::mat xi(T, 1, arma::fill::zeros);
  arma::mat xi_Y(T, 1);
  arma::mat xi_X(T, p);
  double sigma2;
  double IC;
  /*
  arma::mat FE(T, N, arma::fill::zeros) ;
  arma::mat fit(T, N, arma::fill::zeros) ;
  arma::mat FE_add_use(T, N, arma::fill::zeros) ;
  arma::mat FE_inter_use(T, N, arma::fill::zeros) ;
  arma::mat covar_fit(T, N, arma::fill::zeros) ;*/
  arma::mat invXX;

  /* duplicate data */
  arma::mat YY = Y;
  arma::cube XX = X;

  /* grand mean */
  mu_Y = accu(YY) / obs;
  YY = YY - mu_Y;
  if (p > 0) {
    for (int i = 0; i < p; i++) {
      mu_X(i, 0) = accu(XX.slice(i)) / obs;
      XX.slice(i) = XX.slice(i) - mu_X(i, 0);
    }
  }

  /* unit fixed effects */
  if (force == 1 || force == 3) {
    alpha_Y = mean(YY, 0).t();           // colMeans, (N * 1) matrix
    YY = YY - repmat(alpha_Y.t(), T, 1); // (T * N) matrix
    if (p > 0) {
      for (int i = 0; i < p; i++) {
        alpha_X.col(i) = mean(XX.slice(i), 0).t(); // colMeans
        XX.slice(i) = XX.slice(i) - repmat(alpha_X.col(i).t(), T, 1);
      }
    }
  }

  /* time fixed effects  */
  if (force == 2 || force == 3) {
    xi_Y = mean(YY, 1); // rowMeans, (N * 1) matrix
    YY = YY - repmat(xi_Y, 1, N);
    if (p > 0) {
      for (int i = 0; i < p; i++) {
        xi_X.col(i) = mean(XX.slice(i), 1); // rowMeans
        XX.slice(i) = XX.slice(i) - repmat(xi_X.col(i), 1, N);
      }
    }
  }

  /* check if XX has enough variation */
  int p1 = p;
  arma::mat X_invar(p, 1, arma::fill::zeros); // =1 if invar

  int j = 0;
  for (int i = 0; i < p1; i++) {
    if (arma::accu(abs(XX.slice(i))) < 1e-5) {
      XX.shed_slice(i);
      mu_X.shed_row(i);
      alpha_X.shed_col(i);
      xi_X.shed_col(i);
      X_invar(j, 0) = 1;
      i--;
      p1--;
    }
    j++;
  }

  int validX = 1;
  if (p1 == 0) {
    validX = 0;
  } else {
    invXX = XXinv(XX);
  }

  /* Main Algorithm */
  if (p1 == 0) {
    if (r > 0) {
      List pf = panel_factor(YY, r);
      factor = as<arma::mat>(pf["factor"]);
      lambda = as<arma::mat>(pf["lambda"]);
      VNT = as<arma::mat>(pf["VNT"]);
      U = YY - factor * lambda.t();
    } else {
      U = YY;
    }
  } else {
    /* starting value:  the OLS/LSDV estimator */
    if (accu(abs(beta0)) < 1e-10 || r == 0 || b_r != p1) {             //
      beta0 = panel_beta(XX, invXX, YY, arma::zeros<arma::mat>(T, N)); //
    }
    if (r == 0) {
      U = YY;
      beta = beta0;
      for (int k = 0; k < p1; k++) {
        U = U - XX.slice(k) * beta(k, 0);
      }
    } else if (r > 0) {
      // arma::mat invXX =  XXinv(XX) ;  // compute (X'X)^{-1}, outside beta
      // iteration
      List out = beta_iter(XX, invXX, YY, r, tol, beta0, max_iter);
      beta = as<arma::mat>(out["beta"]);
      factor = as<arma::mat>(out["factor"]);
      lambda = as<arma::mat>(out["lambda"]);
      VNT = as<arma::mat>(out["VNT"]);
      U = as<arma::mat>(out["e"]);
      niter = as<int>(out["niter"]);
    }
  }

  /* save fixed effects */
  if (p1 == 0) {
    mu = mu_Y;
    if (force == 1 || force == 3) {
      alpha = alpha_Y;
    }
    if (force == 2 || force == 3) {
      xi = xi_Y;
    }
    /*
    if(r==0){
      List Y_fe_ad ;
      Y_fe_ad = fe_add(alpha_Y, xi_Y, mu_Y, T, N, force) ;
      fit = as<arma::mat>(Y_fe_ad["FE_ad"]) ;
    }
    else if(r>0){
      List E_fe_ad;
      E_fe_ad = fe_add(alpha_Y, xi_Y, mu_Y,T, N, force) ;
      FE_add_use = as<arma::mat>(E_fe_ad["FE_ad"])
      FE_inter_use = factor * lambda.t()
      fit = FE_add_use + FE_inter_use
    }*/
  } else { // with valid covariates

    mu = mu_Y - crossprod(mu_X, beta)(0, 0);

    if (force == 1 || force == 3) {
      alpha = alpha_Y - alpha_X * beta;
    }
    if (force == 2 || force == 3) {
      xi = xi_Y - xi_X * beta;
    }
    /*
    if(r==0){
      for (int i = 0; i < p; i++) {
        covar_fit = covar_fit + XX.slice(i) * beta(i) ;
      }
      List Y_fe_ad ;
      Y_fe_ad = fe_add(alpha, xi, mu, T, N, force) ;
      fit = as<arma::mat>(Y_fe_ad["FE_ad"]) + covar_fit ;
    }
    else if(r>0){
      for (int i = 0; i < p; i++) {
        covar_fit = covar_fit + XX.slice(i) * beta(i) ;
      }
      List E_fe_ad;
      E_fe_ad = fe_add(alpha, xi, mu,T, N, force) ;
      FE_add_use = as<arma::mat>(E_fe_ad["FE_ad"])
      FE_inter_use = factor * lambda.t()
      fit = FE_add_use + FE_inter_use
    }*/
  }

  /* sigma2 and IC */

  // number of estimated parameters
  // force = 0
  double np = r * (N + T) - pow(double(r), 2) + p1 + 1;
  if (force == 1) {
    np = np + (N - 1) - r;
  } else if (force == 2) {
    np = np + (T - 1) - r;
  } else if (force == 3) {
    np = np + (N - 1) + (T - 1) - 2 * r;
  }

  sigma2 = trace(U * U.t()) / (N * T - np);

  IC = log(sigma2) + np * log(double(N * T)) / (N * T);

  // PC criterion in Li 2018
  // mean squared error
  double mse = trace(U * U.t()) / (N * T);

  double m1 = 0;
  if (N < 60) {
    m1 = 60 - N;
  }
  double m2 = 0;
  if (T < 60) {
    m2 = 60 - T;
  }

  double C = (N + m1) * (T + m2) / (N * T);
  double PC = mse + r * sigma2 * C * (N + T) / (N * T) *
                        log(double(N * T) / double(N + T));

  //-------------------------------#
  // Storage
  //-------------------------------#

  List output;

  output["mu"] = mu;
  // output["p1"] = p1 ;

  if (p > 0) {
    // output["beta_valid"] = beta ;
    arma::mat beta_total(p, 1);
    // arma::mat beta_tot(p,1);

    if (p > p1) {
      int j4 = 0;
      for (int i = 0; i < p; i++) {
        if (X_invar(i, 0) == 1) {
          beta_total(i, 0) = arma::datum::nan;
          // beta_tot(i,0) = 0;
        } else {
          beta_total(i, 0) = beta(j4, 0);
          // beta_tot(i,0) = beta(j4,0);
          j4++;
        }
      }
    } else {
      beta_total = beta;
      // beta_tot = beta;
    }
    output["beta"] = beta_total;
    // output["beta_tot"] = beta_tot;
  }

  if (r > 0) {
    output["factor"] = factor;
    output["lambda"] = lambda;
    output["VNT"] = VNT;
    // FE = factor * lambda.t() ;
    // output["FE"] = FE ;
  }
  if ((p1 > 0) && (r > 0)) {
    output["niter"] = niter;
  }
  if (force == 1 || force == 3) {
    output["alpha"] = alpha;
  }
  if (force == 2 || force == 3) {
    output["xi"] = xi;
  }
  output["residuals"] = U;
  output["sigma2"] = sigma2;
  output["IC"] = IC;
  output["PC"] = PC;
  output["validX"] = validX;
  return (output);
}

/* Interactive Fixed Effects: ub */
// [[Rcpp::export]]
List inter_fe_ub(
    const arma::mat& Y, const arma::mat& Y0, const arma::cube& X, const arma::mat& I, const arma::mat& W_in,
    const arma::mat& beta0,
    int r, // r > 0, the outcome has a factor-type fixed effect; r = 0 else
    int force, double tol = 1e-5, int max_iter = 1000) {

  arma::mat W = W_in;
  /* Dimensions */
  int T = Y.n_rows;
  int N = Y.n_cols;
  int p = X.n_slices;
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
  arma::cube XX = X;
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
  if (p1 == 0) {
    if (r > 0) {
      // add fe ; inter fe ; iteration
      List fe_ad_inter =
          fe_ad_inter_iter(YY, Y0, I, W, force, 0, r, 0, 0, tol, max_iter);
      mu = as<double>(fe_ad_inter["mu"]);
      U = as<arma::mat>(fe_ad_inter["e"]);
      fit = as<arma::mat>(fe_ad_inter["fit"]);

      factor = as<arma::mat>(fe_ad_inter["factor"]);
      lambda = as<arma::mat>(fe_ad_inter["lambda"]);
      VNT = as<arma::mat>(fe_ad_inter["VNT"]);
      if (use_weight == 1) {
        burn_in = as<int>(fe_ad_inter["burn_in"]);
      }
      if (force == 1 || force == 3) {
        alpha = as<arma::mat>(fe_ad_inter["alpha"]);
      }
      if (force == 2 || force == 3) {
        xi = as<arma::mat>(fe_ad_inter["xi"]);
      }
      niter = as<int>(fe_ad_inter["niter"]);
    } else {
      if (force == 0) {
        U = YY;
        fit.fill(mu);
      } else {
        // add fe; iteration
        List fe_ad = fe_ad_iter(YY, Y0, I, W, force, tol, max_iter);
        mu = as<double>(fe_ad["mu"]);
        U = as<arma::mat>(fe_ad["e"]);
        fit = as<arma::mat>(fe_ad["fit"]);
        if (force == 1 || force == 3) {
          alpha = as<arma::mat>(fe_ad["alpha"]);
        }
        if (force == 2 || force == 3) {
          xi = as<arma::mat>(fe_ad["xi"]);
        }
        niter = as<int>(fe_ad["niter"]);
      }
    }
  } else {
    /* starting value:  the OLS estimator */
    if (use_weight == 1) {
      invXX = wXXinv(XX, W);
    } else {
      invXX = XXinv(XX); // compute (X'X)^{-1}, outside beta iteration
    }

    if (r == 0) {
      // add fe, covar; iteration
      List fe_ad = fe_ad_covar_iter(XX, invXX, YY, Y0, I, beta0, W, force, tol,
                                    max_iter);
      mu = as<double>(fe_ad["mu"]);
      beta = as<arma::mat>(fe_ad["beta"]);
      U = as<arma::mat>(fe_ad["e"]);
      fit = as<arma::mat>(fe_ad["fit"]);
      if (force == 1 || force == 3) {
        alpha = as<arma::mat>(fe_ad["alpha"]);
      }
      if (force == 2 || force == 3) {
        xi = as<arma::mat>(fe_ad["xi"]);
      }
      niter = as<int>(fe_ad["niter"]);
    } else if (r > 0) {
      // add, covar, interactive, iteration
      List fe_ad_inter_covar = fe_ad_inter_covar_iter(
          XX, invXX, YY, Y0, I, W, beta0, force, 0, r, 0, 0, tol, max_iter);
      mu = as<double>(fe_ad_inter_covar["mu"]);
      beta = as<arma::mat>(fe_ad_inter_covar["beta"]);
      U = as<arma::mat>(fe_ad_inter_covar["e"]);
      fit = as<arma::mat>(fe_ad_inter_covar["fit"]);
      if (use_weight == 1) {
        burn_in = as<int>(fe_ad_inter_covar["burn_in"]);
      }
      factor = as<arma::mat>(fe_ad_inter_covar["factor"]);
      lambda = as<arma::mat>(fe_ad_inter_covar["lambda"]);
      VNT = as<arma::mat>(fe_ad_inter_covar["VNT"]);

      if (force == 1 || force == 3) {
        alpha = as<arma::mat>(fe_ad_inter_covar["alpha"]);
      }
      if (force == 2 || force == 3) {
        xi = as<arma::mat>(fe_ad_inter_covar["xi"]);
      }
      niter = as<int>(fe_ad_inter_covar["niter"]);
    }
  }
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
  return (output);
}
