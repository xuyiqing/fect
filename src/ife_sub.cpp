#include "fect.h"

// sub-functions for unbalanced ife model

/* Obtain additive fe for ub data; assume r=0, without covar */
// [[Rcpp::export]]
List fe_ad_iter (arma::mat Y,
                 arma::mat Y0,
                 arma::mat I,
                 int force,
                 double tolerate) {
  
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  double mu = 0 ;
  double dif = 1.0 ;
  int niter = 0 ;

  arma::mat fit = Y0 ; // initial value
  arma::mat fit_old = Y0 ;
  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;
  double mu_Y = 0 ;
  arma::mat alpha_Y(N, 1, arma::fill::zeros) ;
  arma::mat xi_Y(T, 1, arma::fill::zeros) ;

  arma::mat e(T, N, arma::fill::zeros) ; // residual

  arma::mat YY = Y ;

  List Y_ad ;
  List Y_fe_ad ;

  while (dif > tolerate && niter <= 500) {

    YY =  E_adj (Y, fit, I) ; // e step: expeactation
    
    Y_ad = Y_demean(YY, force) ; 
    mu_Y = as<double>(Y_ad["mu_Y"]) ;
    if (force==1||force==3) {     
      alpha_Y = as<arma::mat>(Y_ad["alpha_Y"]) ;
    }
    if (force==2||force==3) {
      xi_Y = as<arma::mat>(Y_ad["xi_Y"]) ;
    }
    Y_fe_ad = fe_add(alpha_Y, xi_Y, mu_Y, T, N, force) ; // m step: estimate fe
    
    fit = as<arma::mat>(Y_fe_ad["FE_ad"]) ;
    
    dif = arma::norm(fit - fit_old, "fro")/arma::norm(fit_old, "fro") ;
    fit_old = fit ;

    niter = niter + 1 ;
  }

  e = FE_adj(YY - fit, I) ;

  List result;
  mu = as<double>(Y_fe_ad["mu"]) ;
  result["mu"] = mu ;
  result["fit"] = fit ;
  result["niter"] = niter ;
  result["e"] = e ;

  if (force==1||force==3) {
    alpha = as<arma::mat>(Y_fe_ad["alpha"]) ;
    result["alpha"] = alpha ;
  }
  if (force==2||force==3) {
    xi = as<arma::mat>(Y_fe_ad["xi"]) ;
    result["xi"] = xi ;
  }
  return(result) ;
}


/* Obtain additive fe for ub data; assume r=0, with covariates */
// [[Rcpp::export]]
List fe_ad_covar_iter (arma::cube XX,
                       arma::mat xxinv,
                       arma::mat Y,
                       arma::mat Y0,
                       arma::mat I,
                       arma::mat beta0,
                       int force,
                       double tolerate) {
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  int p = XX.n_slices ;
  double dif = 1.0 ;
  int niter = 0 ;
  
  arma::mat beta(p, 1, arma::fill::zeros) ;
  if (beta0.n_rows == p) {
    beta = beta0 ;
  }

  double mu = 0 ;

  arma::mat fit(T, N, arma::fill::zeros) ;
  arma::mat fit_old(T, N, arma::fill::ones) ;

  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;

  arma::mat YY = Y ;

  arma::mat U(T, N, arma::fill::zeros) ; 
  arma::mat e(T, N, arma::fill::zeros) ; // residual

  List ife_inner ;

  arma::mat covar_fit(T, N, arma::fill::zeros) ;
  for (int i = 0; i < p; i++) {
    covar_fit = covar_fit + XX.slice(i) * beta(i) ;
  }

  fit = Y0 ;
  fit_old = fit ;

  arma::mat FE = fit - covar_fit ; // initial fixed effects

  while (dif > tolerate && niter <= 500) {

    YY =  E_adj (Y, fit, I) ; // e-step: expectation
    
    // m1: estimate beta
    beta = panel_beta(XX, xxinv, YY, FE) ;
    covar_fit.zeros() ;
    for (int i = 0; i < p; i++) {
      covar_fit = covar_fit + XX.slice(i) * beta(i) ;
    }

    // m2: estimate interactive fe, additive fe, and mu
    U = E_adj (YY - covar_fit, FE, I) ;
    ife_inner = ife(U, force, 0, 0, 0, 0) ;
    FE = as<arma::mat>(ife_inner["FE"]) ; // overall fe 

    fit = covar_fit + FE ;

    dif = arma::norm(fit - fit_old, "fro")/arma::norm(fit_old, "fro") ;
    fit_old = fit ;

    niter = niter + 1 ;
  }

  e = FE_adj(YY - fit, I) ;

  List result ;
  mu = as<double>(ife_inner["mu"]) ;
  result["mu"] = mu ;
  result["fit"] = fit ;
  result["niter"] = niter ;
  result["e"] = e ;
  if (p>0) {
    result["beta"] = beta ;
  }
  if (force==1||force==3) {
    alpha = as<arma::mat>(ife_inner["alpha"]) ;
    result["alpha"] = alpha ;
  }
  if (force==2||force==3) {
    xi = as<arma::mat>(ife_inner["xi"]) ;
    result["xi"] = xi ;
  }
  return(result) ;
}

/* Obtain additive fe for ub data; assume r>0 but p=0*/
// [[Rcpp::export]]
List fe_ad_inter_iter (arma::mat Y,
                       arma::mat Y0,
                       arma::mat I,
                       int force,
                       int mc, // whether pac or mc method
                       int r,
                       int hard, 
                       double lambda,
                       double tolerate
                       ) {
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  double mu = 0 ;
  double dif = 1.0 ;
  int niter = 0 ;
  int validF = 1 ; // whether has a factor structure

  arma::mat VNT(r, r) ;
  arma::mat FE_inter_use(T, N, arma::fill::zeros) ;
  arma::mat fit(T, N, arma::fill::zeros) ;
  arma::mat fit_old(T, N, arma::fill::ones) ;
  arma::mat e(T, N, arma::fill::zeros) ; // residual
  arma::mat U(T, N, arma::fill::zeros) ;
  
  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;

  arma::mat F(T, r, arma::fill::zeros) ;
  arma::mat L(N, r, arma::fill::zeros) ;


  arma::mat YY = Y ;

  List pf ;
  List ife_inner ;

  // initial value for ife

  //if (hard == 0) {
  //  U = FE_adj(Y - Y0, I) ;
  //  if (mc == 0) {
  //    pf = panel_factor(U, r)  ;
  //    F = as<arma::mat>(pf["factor"]) ;
  //    L = as<arma::mat>(pf["lambda"]) ;
  //    FE_inter_use = F * L.t() ; // interactive fe
  //  }
  //  else {
  //    FE_inter_use = panel_FE(U, lambda, hard) ;
  //  }
  //  fit = Y0 + FE_inter_use ;
  //  fit_old = fit ;
  //} else {
    fit = Y0 ;
    fit_old = fit ;
  //}
  
  while (dif > tolerate && niter <= 500) {
    
    YY = E_adj(Y, fit, I) ; // e-step: expectation
    
    if (mc == 0) {
      ife_inner = ife(YY, force, 0, r, 0, 0) ;
    }
    else {
      ife_inner = ife(YY, force, 1, 1, hard, lambda) ;
    }
    fit = as<arma::mat>(ife_inner["FE"]) ; // new overall fe

    dif = arma::norm(fit - fit_old, "fro")/arma::norm(fit_old, "fro") ;
    fit_old = fit ;

    niter = niter + 1 ;
  }
  e = FE_adj(YY - fit, I) ;
  FE_inter_use = as<arma::mat>(ife_inner["FE_inter_use"]) ;
  if (arma::accu(abs(FE_inter_use)) < 1e-10) {
    validF = 0 ;
  }

  List result;
  mu = as<double>(ife_inner["mu"]) ;
  result["mu"] = mu ;
  result["niter"] = niter ;
  result["fit"] = fit ;
  result["e"] = e ;
  result["validF"] = validF ;
  if (force==1||force==3) {
    alpha = as<arma::mat>(ife_inner["alpha"]) ;
    result["alpha"] = alpha ;
  }
  if (force==2||force==3) {
    xi = as<arma::mat>(ife_inner["xi"]) ;
    result["xi"] = xi ;
  }
  if (mc == 0) {
    result["lambda"] = as<arma::mat>(ife_inner["lambda"]) ;
    result["factor"] = as<arma::mat>(ife_inner["factor"]) ;
    result["VNT"] = as<arma::mat>(ife_inner["VNT"]) ;
  }
  return(result) ;
}

/* Obtain additive fe for ub data; assume r>0 p>0*/
// [[Rcpp::export]]
List fe_ad_inter_covar_iter (arma::cube XX,
                             arma::mat xxinv,
                             arma::mat Y,
                             arma::mat Y0,
                             arma::mat I,
                             arma::mat beta0, 
                             int force,
                             int mc, // whether pca or mc method
                             int r,
                             int hard, 
                             double lambda,
                             double tolerate
                             ) {
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  int p = XX.n_slices ;
  double dif = 1.0 ;
  int niter = 0 ;
  int validF = 1 ;

  arma::mat beta(p, 1, arma::fill::zeros) ;
  if (beta0.n_rows == p) {
    beta = beta0 ;
  }

  double mu = 0 ;
  arma::mat VNT(r, r) ;
  arma::mat FE_inter_use(T, N, arma::fill::zeros) ; // ife

  arma::mat fit(T, N, arma::fill::zeros) ;
  arma::mat fit_old(T, N, arma::fill::ones) ;
  arma::mat U(T, N, arma::fill::zeros) ;

  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;

  arma::mat YY = Y ;

  arma::mat e(T, N, arma::fill::zeros) ; // residual

  List ife_inner ;
  List pf ;

  arma::mat F ;
  arma::mat L ;

  // initial value for ife
  
  //if (hard == 0) {
  //  U = FE_adj(Y - Y0, I) ;
  //  if (mc == 0) {
  //    pf = panel_factor(U, r)  ;
  //    F = as<arma::mat>(pf["factor"]) ;
  //    L = as<arma::mat>(pf["lambda"]) ;
  //    FE_inter_use = F * L.t() ; // interactive fe
  //  }
  //  else {
  //    FE_inter_use = panel_FE(U, lambda, hard) ;
  //  }

  //  fit = Y0 + FE_inter_use ;
  //  fit_old = fit ;
  //} else {
    fit = Y0 ;
    fit_old = fit ;
  //}

  arma::mat covar_fit(T, N, arma::fill::zeros) ;
  for (int i = 0; i < p; i++) {
    covar_fit = covar_fit + XX.slice(i) * beta(i) ;
  }
  arma::mat FE = fit - covar_fit ;
  
  while (dif > tolerate && niter <= 500) {
    
    YY =  E_adj (Y, fit, I) ; // e-step: expectation
    
    // m1: estimate beta
    beta = panel_beta(XX, xxinv, YY, FE) ;
    
    covar_fit.zeros() ;
    for (int i = 0; i < p; i++) {
      covar_fit = covar_fit + XX.slice(i) * beta(i) ;
    }
    
    // m2: estimate interactive fe, additive fe, and mu
    U = E_adj (YY - covar_fit, FE, I) ;
    ife_inner = ife(U, force, mc, r, hard, lambda) ;
    FE = as<arma::mat>(ife_inner["FE"]) ;

    fit = covar_fit + FE ; // overall fe */

    dif = arma::norm(fit - fit_old, "fro")/arma::norm(fit_old, "fro") ;
    fit_old = fit ;

    niter = niter + 1 ;
  }
  e = FE_adj(Y - fit, I) ;

  FE_inter_use = as<arma::mat>(ife_inner["FE_inter_use"]) ;
  if (arma::accu(abs(FE_inter_use)) < 1e-10) {
    validF = 0 ;
  }

  List result;
  mu = as<double>(ife_inner["mu"]) ;
  result["mu"] = mu ;
  result["niter"] = niter ;
  result["e"] = e ;
  result["beta"] = beta ;
  result["fit"] = fit ;
  result["validF"] = validF ;

  if (force==1||force==3) {
    alpha = as<arma::mat>(ife_inner["alpha"]) ;
    result["alpha"] = alpha ;
  }
  if (force==2||force==3) {
    xi = as<arma::mat>(ife_inner["xi"]) ;
    result["xi"] = xi ;
  }
  if (mc == 0) {
    L = as<arma::mat>(ife_inner["lambda"]) ;
    F = as<arma::mat>(ife_inner["factor"]) ;
    VNT = as<arma::mat>(ife_inner["VNT"]) ;
    result["lambda"] = L ;
    result["factor"] = F ;
    result["VNT"] = VNT ;
  }
  return(result) ;
}

/* Main iteration for beta */
// [[Rcpp::export]]
List beta_iter (arma::cube X,
                arma::mat xxinv,
                arma::mat Y,
                int r,
                double tolerate,
                arma::mat beta0) {

  /* beta.new: computed beta under iteration with error precision=tolerate
     factor: estimated factor
     lambda: estimated loadings
     V: the eigenvalues matrix
     e: estimated residuals
     niter: number of interations to achieve convergence */
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  int p = X.n_slices ;
  int b_r = beta0.n_rows ; 
  double beta_norm = 1.0 ;
  arma::mat beta(p, 1, arma::fill::zeros) ;
  if (b_r == p) {
      beta = beta0 ;
  } // beta should have the same dimension as X, if not it will be reset to 0
  arma::mat beta_old = beta ;
  arma::mat VNT(r, r, arma::fill::zeros) ;
  arma::mat FE(T, N, arma::fill::zeros) ;

  /* starting value */
  arma::mat U = Y ;
  for (int k = 0; k < p; k++) {
    U = U - X.slice(k) * beta(k) ;
  }
  List pf = panel_factor(U, r)  ;
  arma::mat F = as<arma::mat>(pf["factor"]) ;
  arma::mat L = as<arma::mat>(pf["lambda"]) ;
 
  /* Loop */
  int niter = 0 ;
  while ((beta_norm > tolerate) && (niter < 500)) {
    niter++ ; 
    FE = F * L.t() ;
    beta = panel_beta(X, xxinv, Y, FE) ;
    beta_norm = arma::norm(beta - beta_old, "fro") ; 
    beta_old = beta ;
    U = Y ;
    for (int k = 0; k < p; k++) {
      U = U - X.slice(k) * beta(k) ;
    }
    pf = panel_factor(U, r)  ;
    F = as<arma::mat>(pf["factor"]) ;
    L = as<arma::mat>(pf["lambda"]) ; 
  }
  VNT = as<arma::mat>(pf["VNT"]) ; 
  arma::mat e = U - F * L.t() ;

  /* Storage */
  List result ;
  result["niter"] = niter ;
  result["beta"] = beta ;
  result["e"] = e ; 
  result["lambda"] = L ;
  result["factor"] = F ;
  result["VNT"] = VNT ;
  return(result)  ;
}
