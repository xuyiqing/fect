#include "fect.h"

// matrix completion

/* Interactive Fixed Effects: matrix completion */
// [[Rcpp::export]]
List inter_fe_mc (arma::mat Y,
                  arma::mat Y0,
                  arma::cube X,
                  arma::mat I,
                  arma::mat beta0,
                  int r, // r > 0, the outcome has a factor-type fixed effect; r = 0 else
                  double lambda,
                  int force,
                  double tol = 1e-5
                  ) {
  
  /* Dimensions */
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  int p = X.n_slices ;
  double obs = accu(I) ;
  int niter = 0 ;
  int validF = 1 ;
  arma::mat beta ; 
  arma::mat U ;
  double mu = 0 ;
  double mu_Y = 0 ;

  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;
  arma::mat fit(T, N, arma::fill::zeros) ;
  // double sigma2 = 0;


  arma::mat invXX ;

  /* duplicate data */
  arma::mat YY = Y;
  arma::cube XX = X;

  
  /* check if XX has enough variation */
  int p1 = p; 
  arma::mat X_invar(p, 1, arma::fill::zeros); // =1 if invar

  int j = 0;

  for(int i=0; i<p1; i++){
    if (arma::accu(abs(XX.slice(i))) < 1e-5) {
      XX.shed_slice(i);
      X_invar(j,0) = 1;
      i--;
      p1--;
    }
    j++;
  }

  int validX = 1 ;
  if(p1==0){
    validX = 0 ;
    if (force == 0 && r == 0) { // no covariate and force == 0 and r == 0 
      mu_Y = accu(YY)/obs ;
      mu = mu_Y ;
      YY = FE_adj(YY - mu_Y, I) ;
    }
  }

  /* Main Algorithm */ 
  if (p1 == 0) {
    if (r > 0) {
      // add fe ; inter fe ; iteration

      // soft impute as starting value
      List fe_ad_inter = fe_ad_inter_iter(YY, Y0, I, force, 1, 1, 0, lambda, tol) ;
      
      // fit = as<arma::mat>(fe_ad_inter["fit"]) ;
      // hard impute 
      // fe_ad_inter = fe_ad_inter_iter(YY, fit, I, force, 1, 1, 1, lambda, tol) ;
      
      fit = as<arma::mat>(fe_ad_inter["fit"]) ;
      mu = as<double>(fe_ad_inter["mu"]) ;
      U = as<arma::mat>(fe_ad_inter["e"]) ;
      if (force==1||force==3) {
        alpha = as<arma::mat>(fe_ad_inter["alpha"]) ;
      }
      if (force==2||force==3) {
        xi = as<arma::mat>(fe_ad_inter["xi"]) ;
      }
      niter = as<int>(fe_ad_inter["niter"]) ;
      validF = as<int>(fe_ad_inter["validF"]) ;
    } 
    else {
      if (force==0) {
        U = YY ;
        fit.fill(mu) ;
        validF = 0 ;
      } else {
        // add fe; iteration
        List fe_ad = fe_ad_iter(YY, Y0, I, force, tol) ;
        mu = as<double>(fe_ad["mu"]) ;
        U = as<arma::mat>(fe_ad["e"]) ;
        fit = as<arma::mat>(fe_ad["fit"]) ;
        if (force==1||force==3) {
          alpha = as<arma::mat>(fe_ad["alpha"]) ;
        }
        if (force==2||force==3) {
          xi = as<arma::mat>(fe_ad["xi"]) ;
        }
        niter = as<int>(fe_ad["niter"]) ;
        validF = 0 ;
      }
    } 
  } 
  else {
    /* starting value:  the OLS estimator */
    invXX = XXinv(XX) ; // compute (X'X)^{-1}, outside beta iteration 
    if (r==0) {
      // add fe, covar; iteration
      List fe_ad = fe_ad_covar_iter(XX, invXX,
                                    YY, Y0, I, beta0, force, tol) ;
      mu = as<double>(fe_ad["mu"]) ;
      beta = as<arma::mat>(fe_ad["beta"]) ;
      U = as<arma::mat>(fe_ad["e"]) ;
      fit = as<arma::mat>(fe_ad["fit"]) ;
      if (force==1||force==3) {
        alpha = as<arma::mat>(fe_ad["alpha"]) ;
      }
      if (force==2||force==3) {
        xi = as<arma::mat>(fe_ad["xi"]) ;
      }
      niter = as<int>(fe_ad["niter"]) ;
      validF = 0 ;
    } 
    else if (r > 0) {       
      // add, covar, interactive, iteration
      // soft impute as starting value
      List fe_ad_inter_covar = fe_ad_inter_covar_iter(XX, invXX,
               YY, Y0, I, beta0, 
               force, 1, 1, 0, lambda, tol) ;
      beta = as<arma::mat>(fe_ad_inter_covar["beta"]) ;
      fit = as<arma::mat>(fe_ad_inter_covar["fit"]) ;

      // hard impute 
      // fe_ad_inter_covar = fe_ad_inter_covar_iter(XX, invXX,
      //           YY, fit, I, beta, 
      //           force, 1, 1, 1, lambda, tol) ;
      
      // beta = as<arma::mat>(fe_ad_inter_covar["beta"]) ;
      // fit = as<arma::mat>(fe_ad_inter_covar["fit"]) ;
      mu = as<double>(fe_ad_inter_covar["mu"]) ;
      U = as<arma::mat>(fe_ad_inter_covar["e"]) ;
      if (force==1||force==3) {
        alpha = as<arma::mat>(fe_ad_inter_covar["alpha"]) ;
      }
      if (force==2||force==3) {
        xi = as<arma::mat>(fe_ad_inter_covar["xi"]) ;
      }
      niter = as<int>(fe_ad_inter_covar["niter"]) ;
      validF = as<int>(fe_ad_inter_covar["validF"]) ;
    }
  } 
    
  //-------------------------------#
  // Storage
  //-------------------------------# 

  List output ;

  if(p>0) {
    arma::mat beta_total(p,1);

    if(p>p1) {
      int j4= 0;
      for(int i=0; i<p; i++) {
        if(X_invar(i,0)==1) {
          beta_total(i,0) = arma::datum::nan;
        }
        else {
          beta_total(i,0) = beta(j4,0);
          j4++;
        }
      }
    }
    else {
      beta_total = beta;
    }
    output["beta"] = beta_total;
  }

  output["mu"] = mu ;   
  output["fit"] = fit ; 
  output["validF"] = validF ; 

  // if ( !(force == 0 && r == 0 && p1 == 0) ) {
    output["niter"] = niter ;
  // }
  if (force ==1 || force == 3) {
    output["alpha"] = alpha ;
  }
  if (force ==2 || force == 3) {
    output["xi"] = xi ;
  }
  output["residuals"] = U ;
  output["validX"] = validX;
  return(output);
 
}
