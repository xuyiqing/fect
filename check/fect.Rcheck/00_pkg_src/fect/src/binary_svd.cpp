#include "fect.h"

// core binary function

            /* ---------------------- */
/* ----------- binary ife using SVD ------------- */
            /* ---------------------- */

/* inter fe model for binary outcome */
// [[Rcpp::export]]
List inter_fe_d (arma::mat Y,
                 arma::mat Y_fit0, // initial fit using ols
                 arma::mat FE0, // initial fixed effects 
                 arma::cube X,
                 int r,
                 int force,
                 int mniter = 5000,
                 double w = 1.0,
                 double tol = 1e-5) {
  
  /* Dimensions */
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  int p = X.n_slices ;
  int niter = 0 ;
  int validX = 1 ;
  double dif = 1.0 ;

  double allh_old = 0 ;
  double allh ;
  double llh ;
  double IC ;

  double mu = 0 ;
  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat lambda(N, r, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ; 
  arma::mat factor(T, r, arma::fill::zeros) ;
  arma::mat VNT(r, r, arma::fill::zeros) ;
  
  arma::mat YY = Y ;
  arma::cube XX = X ;
  arma::mat invXX(p, p, arma::fill::zeros) ;
  arma::mat beta(p, 1, arma::fill::zeros) ;
  invXX =  XXinv(XX) ;

  arma::mat M ;
  arma::mat Y_fit_x(T, N, arma::fill::zeros) ;
  arma::mat FE(T, N, arma::fill::zeros) ;
  arma::mat Y_star(T, N, arma::fill::zeros) ;
  arma::mat res(T, N, arma::fill::zeros) ;
  arma::mat V_iter ;

  // initialize
  arma::mat Y_fit = Y_fit0 ;
  arma::mat Y_fit_old = Y_fit ;
  FE = FE0 ;

  List ife_inner ;

  double g_h = 0 ;
  double s_mo = 0 ;


  while ((dif>tol) && (niter < mniter)) {
    niter++ ;
    
    M = M_gen(Y_fit, Y) ;
    Y_star = Y_fit + M ; // expectation
    
    // update beta
    if (niter > 1) {
      if (p > 0) {
        beta = panel_beta(XX, invXX, Y_star, FE) ;
        Y_fit_x.zeros() ;
        for (int i = 0; i < p; i++) {
          Y_fit_x = Y_fit_x + XX.slice(i) * beta(i) ;
        }
      }
    }

    // update fixed effects
    res = Y_star - Y_fit_x ;  
    ife_inner = ife(res, force, 0, r, 0, 0) ;
    FE = as<arma::mat>(ife_inner["FE"]) ; // interactive fe

    // latent outcome and residuals
    Y_fit = Y_fit_x + FE ;
    res = Y_star - Y_fit ;

    // monotonically overrelaxed
    V_iter = V(Y, Y_fit_old) ; // impute sigma
    g_h = gamma_hat(res, V_iter) ;
    s_mo = S(1.0, g_h, w) ;

    // final adjustment
    Y_fit = Y_fit * pow(s_mo,0.5) ;
    Y_fit_old = Y_fit ;
    FE = FE * pow(s_mo, 0.5) ;

    llh = loglh(Y_fit,Y) ; 
    allh = llh/(N*T) ;
    dif = std::abs(allh - allh_old) ;
    allh_old = allh ;
  }
  // BIC
  IC = (r*(N+T)-pow(double(r),2)+p)*log(double(N*T))/(N*T)-2*llh/(N*T) ;

  //-----------------##
  // storage
  //-----------------##
  mu = as<double>(ife_inner["mu"]) ;
  mu = mu * pow(s_mo,0.5) ;
  
  List output ;
  output["niter"] = niter ; 
  output["fit"] = Y_fit ;
  output["mu"] = mu;
  output["loglikelihood"] = allh ;
  
  if ( validX == 1 ) {
    output["beta"] = beta * pow(s_mo,0.5) ;
  }  
  output["validX"]=validX;
  output["IC"]=IC;
  
  if (r > 0) {
    factor=as<arma::mat>(ife_inner["factor"]);
    lambda=as<arma::mat>(ife_inner["lambda"]);
    VNT=as<arma::mat>(ife_inner["VNT"]);
    lambda = lambda * pow(s_mo,0.5);
    VNT = VNT * pow(s_mo,0.5);
    output["factor"] = factor ;
    output["lambda"] = lambda ;
    output["VNT"] = VNT ;
  }
  if( force == 1 || force == 3 ){
    alpha = as<arma::mat>(ife_inner["alpha"]) ;
    alpha = alpha * pow(s_mo,0.5) ;
    output["alpha"] = alpha ;
  }
  if( force == 2 || force == 3 ){
    xi = as<arma::mat>(ife_inner["xi"]) ;
    xi = xi * pow(s_mo,0.5) ;
    output["xi"] = xi ;
  }
  return(output);
}

/* inter fe model for binary outcome: unbalanced panel */
// [[Rcpp::export]]
List inter_fe_d_ub (arma::mat Y, 
                    arma::mat Y_fit0, // initial fit using ols
                    arma::mat FE0, // initial fixed effects
                    arma::cube X,
                    arma::mat I,
                    int r,
                    int force, 
                    int mniter = 5000,
                    double w = 1.0,
                    double tol = 1e-5) {
  /* Dimensions */
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  int p = X.n_slices ;
  int obs = arma::accu(I) ;
  int niter = 0 ;
  int validX = 1 ;
  double dif = 1.0 ;
  double mu = 0 ;
  double allh_old = 0;
  double allh;
  double llh;
  double IC;
  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat lambda(N, r, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros);
  arma::mat factor(T, r, arma::fill::zeros) ;

  arma::mat VNT(r, r, arma::fill::zeros) ;
  
  arma::mat YY = Y ;
  arma::cube XX = X ;
  arma::mat invXX(p, p, arma::fill::zeros) ;
  arma::mat beta(p, 1, arma::fill::zeros) ;


  arma::mat M ;
  arma::mat Y_fit(T, N, arma::fill::zeros) ;
  arma::mat Y_fit_old(T, N, arma::fill::zeros) ;
  arma::mat Y_fit_x(T, N, arma::fill::zeros) ;
  arma::mat FE(T, N, arma::fill::zeros) ;
  arma::mat Y_star(T, N, arma::fill::zeros) ;
  arma::mat res(T, N, arma::fill::zeros) ;
  arma::mat X_sub(T, N, arma::fill::zeros) ;
  arma::mat V_iter ;

  // initialize
  Y_fit = Y_fit0 ;
  Y_fit_old = Y_fit ;
  FE = FE0 ;

  List ife_inner ;

  double g_h = 0 ;
  double s_mo = 0 ;

  if (p > 0) {
    for (int i = 0; i < p; i++) {
      X_sub = XX.slice(i) ;
      XX.slice(i) = FE_adj(X_sub, I) ;
    }
    invXX =  XXinv(XX) ;
  } else {
    validX = 0 ;
  }

  while ((dif>tol) && (niter < mniter)) {
    niter++ ;
    M = M_gen_ub(Y_fit, Y, I) ;
    Y_star = Y_fit + M ; // expectation
    
    // update beta 
    if (p > 0) {
      beta = panel_beta(XX, invXX, Y_star, FE) ;
      Y_fit_x.zeros() ;
      for (int i = 0; i < p; i++) {
        Y_fit_x = Y_fit_x + XX.slice(i) * beta(i) ;
      }
    }

    // update fixed effects
    res = E_adj(Y_star - Y_fit_x, FE, I) ;
    ife_inner = ife(res, force, 0, r, 0, 0) ;
    FE = as<arma::mat>(ife_inner["FE"]) ; // interactive fe
    
    Y_fit = Y_fit_x + FE ;
    res = FE_adj(Y_star - Y_fit, I) ;

    // monotonically overrelaxed
    V_iter = V_ub(Y, Y_fit_old, I) ; // impute sigma
    g_h = gamma_hat_ub(res, V_iter, I) ;
    s_mo = S(1.0, g_h, w) ;

    // final adjustment
    Y_fit = Y_fit * pow(s_mo,0.5) ;
    Y_fit_old = Y_fit ;
    FE = FE * pow(s_mo, 0.5) ;

    llh = loglh_ub(Y_fit, Y, I) ; 
    allh = llh/obs ;
    dif = std::abs(allh - allh_old);
    allh_old = allh;
  }
  // BIC
  IC = (r*(N+T)-pow(double(r),2)+p)*log(obs)/(obs)-2*llh/(obs) ;

  //-----------------##
  // storage
  //-----------------##
  mu = as<double>(ife_inner["mu"]) ;
  mu = mu * pow(s_mo,0.5) ;
  
  List output ;
  output["niter"] = niter ; 
  output["mu"] = mu;
  output["loglikelihood"] = allh ;
  
  if ( validX == 1 ) {
    beta = beta * pow(s_mo,0.5) ;
    Y_fit = FE ; 
    for (int i = 0; i < p; i++) {
      Y_fit = Y_fit + X.slice(i) * beta(i) ;
    }
    output["beta"] = beta ;
  } 
  output["fit"] = Y_fit ;
  output["validX"]=validX;
  output["IC"]=IC;
  
  if (r > 0) {
    factor=as<arma::mat>(ife_inner["factor"]);
    lambda=as<arma::mat>(ife_inner["lambda"]);
    VNT=as<arma::mat>(ife_inner["VNT"]);
    lambda = lambda * pow(s_mo,0.5);
    VNT = VNT * pow(s_mo,0.5);
    output["factor"] = factor ;
    output["lambda"] = lambda ;
    output["VNT"] = VNT ;
  }
  if( force == 1 || force == 3 ){
    alpha = as<arma::mat>(ife_inner["alpha"]) ;
    alpha = alpha * pow(s_mo,0.5) ;
    output["alpha"] = alpha ;
  }
  if( force == 2 || force == 3 ){
    xi = as<arma::mat>(ife_inner["xi"]) ;
    xi = xi * pow(s_mo,0.5) ;
    output["xi"] = xi ;
  }
  return(output);
}
