#include "fect.h"

// sub-functions for binary model

/* Obtain factors and loadings by QR decomposition */
// [[Rcpp::export]]
List qr_factor (const arma::mat& F, const arma::mat& L) {
  
  int T = F.n_rows ;
  int N = L.n_rows ;
  int r = F.n_cols ;

  arma::mat factor(T, r, arma::fill::zeros) ;
  arma::mat lambda(N, r, arma::fill::zeros) ;
  arma::mat FE (T, N, arma::fill::zeros) ;

  arma::mat Q1 ; 
  arma::mat R1 ;
  arma::qr(Q1, R1, F) ;
  arma::mat Qf = Q1.cols(0, r-1) ;
  arma::mat Rf = R1.rows(0, r-1) ;

  // arma::mat RL = Rf * L.t() ;
  // arma::mat RL = L * Rf.t() ;

  arma::mat Q2 ;
  arma::mat R2 ;
  arma::qr(Q2, R2, Rf * L.t()) ;
  // arma::mat Ql = Q2.cols(0, r-1) ;
  // arma::mat Rl = R2.rows(0, r-1) ;

  factor = sqrt(double(T)) * Qf * Q2 ;
  // factor = Qf * Rfl.t() ;
  lambda = R2.t() / sqrt(double(T)) ;
  // factor = sqrt(double(T)) * Qf ;
  // lambda = L * Rf.t() / sqrt(double(T)) ;

  FE = factor * lambda.t() ;
  List result ;
  result["lambda"] = lambda ;
  result["factor"] = factor ;
  result["FE"] = FE ;
  return(result) ;
  
}

/* I a T*1 vector */
// [[Rcpp::export]]
arma::mat IND (const arma::mat& I) {
  int T = I.n_rows ;
  int n = (int) arma::accu(I) ;
  arma::mat ind(n, 1, arma::fill::zeros) ;
  int j = 0 ;
  
  for (int i = 0; i < T; i++) {
    if (I(i,0) == 1) {
      ind(j,0) = i ;
      j = j + 1 ;
    }
  }
  return(ind) ;
}

/* for a given unit or period */
// [[Rcpp::export]]
List subfe (const arma::mat& Y, 
            const arma::mat& X, 
            const arma::mat& I,
            int intercept) {

  int r = X.n_cols ;

  arma::mat ind = IND(I) ;
  arma::mat subY = subr(Y, ind) ; 
  arma::mat subX = subr(X, ind) ;  
  arma::mat beta(r, 1, arma::fill::zeros) ;
  double mu = 0 ;
  int N = subY.n_rows ;
  arma::mat en(N, 1, arma::fill::ones) ;

  arma::mat workingbeta ;
  if (intercept == 0) {
    beta = arma::pinv(subX) * subY ;
  } else {
    subX.insert_cols(0, en) ;
    workingbeta = arma::pinv(subX) * subY ;
    mu = workingbeta(0, 0) ;
    workingbeta.shed_row(0) ;
    beta = workingbeta ;
  }

  List result ;
  result["beta"] = beta.t() ;
  result["mu"] = mu ;
  return(result) ;
}

/* for all units */
// [[Rcpp::export]]
List l_ub (const arma::mat& Y, 
           const arma::mat& F, 
           const arma::mat& I,
           int r,
           int force) {

  // int T = Y.n_rows ;
  int N = Y.n_cols ;
  int intercept = 0 ;

  arma::mat lambda(N, r, arma::fill::zeros) ;
  arma::mat alpha(N, 1, arma::fill::zeros) ;


  List sublambda ;
  if (r > 0) {
    for (int i = 0; i < N; i++) {
      if (r > 0) {
        if (force == 0) {
          sublambda = subfe(Y.col(i), F, I.col(i), 0) ;
          lambda.row(i) = as<arma::mat>(sublambda["beta"]) ;
        }
        else if (force == 1) {
          // alpha_1 = 0 
          intercept = 1 ;
          if (i == 0) {
            intercept = 0 ;
          }
          sublambda = subfe(Y.col(i), F, I.col(i), intercept) ;
          lambda.row(i) = as<arma::mat>(sublambda["beta"]) ;
          if (i == 0) {
            alpha(i, 0) = 0 ;
          } else {
            alpha(i, 0) = as<double>(sublambda["mu"]) ;
          }

        }
        else if (force == 2) {
          // lambda_1 = 0 
          if (i == 0) {
            lambda.row(i) = arma::zeros(1, r) ;
          } else {
            sublambda = subfe(Y.col(i), F, I.col(i), 0) ;
            lambda.row(i) = as<arma::mat>(sublambda["beta"]) ;
          }
        }
        else if (force == 3) {
          // alpha_1 = 0, lambda_1 = 0
          if (i == 0) {
            lambda.row(i) = arma::zeros(1, r) ;
            alpha(i, 0) = 0 ;
          } else {
            sublambda = subfe(Y.col(i), F, I.col(i), 1) ;
            lambda.row(i) = as<arma::mat>(sublambda["beta"]) ;
            alpha(i, 0) = as<double>(sublambda["mu"]) ;
          } 
        }
      } 
    }
  } 
  else {
    if (force == 1 || force == 3) {
      alpha(0, 0) = 0 ;
      for (int i = 1; i < N; i++) {
        alpha(i,0) = accu(Y.col(i))/accu(I.col(i)) ;
      }
    }
  }

  List result ;
  if (r > 0) {
    result["L"] = lambda ;
  }
  if (force == 1 || force == 3) {
    result["alpha"] = alpha ;
  }
  return(result) ;
}


/* for all periods */
// [[Rcpp::export]]
List f_ub (const arma::mat& Y, 
           const arma::mat& L, 
           const arma::mat& I,
           int r,
           int force) {

  int T = Y.n_rows ;
  // int N = Y.n_cols ;
  int intercept = 0 ;

  arma::mat factor(T, r, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;

  List subfactor ;
  if (r > 0) {
    for (int i = 0; i < T; i++) {
      if (r > 0) {
        if (force == 0) {
          subfactor = subfe((Y.row(i)).t(), L, (I.row(i)).t(), 0) ;
          factor.row(i) = as<arma::mat>(subfactor["beta"]) ;
        }
        else if (force == 1) {
          // f_1 = 0
          if (i == 0) {
            factor.row(i) = arma::zeros(1, r) ;
          } else {
            subfactor = subfe((Y.row(i)).t(), L, (I.row(i)).t(), 0) ;
            factor.row(i) = as<arma::mat>(subfactor["beta"]) ;
          }
        }
        else if (force == 2) {
          // xi_1 = 0 
          intercept = 1 ;
          if (i == 0) {
            intercept = 0 ;
          }
          subfactor = subfe((Y.row(i)).t(), L, (I.row(i)).t(), intercept) ;
          factor.row(i) = as<arma::mat>(subfactor["beta"]) ;
          if (i == 0) {
            xi(i, 0) = 0 ;
          } else {
            xi(i, 0) = as<double>(subfactor["mu"]) ;
          }
        }
        else if (force == 3) {
          // xi_1 = 0, f_1 = 0
          if (i == 0) {
            factor.row(i) = arma::zeros(1, r) ;
            xi(i, 0) = 0 ;
          } else {
            subfactor = subfe((Y.row(i)).t(), L, (I.row(i)).t(), 1) ;
            factor.row(i) = as<arma::mat>(subfactor["beta"]) ;
            xi(i, 0) = as<double>(subfactor["mu"]) ;
          } 
        }
      } 
    }
  } 
  else {
    if (force == 2 || force == 3) {
      xi(0, 0) = 0 ;
      for (int i = 1; i < T; i++) {
        xi(i,0) = accu(Y.row(i))/accu(I.row(i)) ;
      }
    }
  }

  List result ;
  if (r > 0) {
    result["F"] = factor ;
  }
  if (force == 2 || force == 3) {
    result["xi"] = xi ;
  }
  return(result) ;
}

/* factor analysis: mu add ife*/
// [[Rcpp::export]]
List fe (const arma::mat& E,
         const arma::mat& F_old_in,
         const arma::mat& xi_old,
         int force,
         int r
        ) {
  
  arma::mat F_old = F_old_in;
  int T = E.n_rows ; 
  int N = E.n_cols ;

  arma::mat et(T, 1, arma::fill::ones) ;
  arma::mat en(N, 1, arma::fill::ones) ;
  
  // complete fixed effects 
  arma::mat u0(1, 1, arma::fill::zeros) ;
  arma::mat f0(1, r, arma::fill::zeros) ;
 
  arma::mat F(T, r, arma::fill::zeros) ;
  arma::mat L(N, r, arma::fill::zeros) ;
  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;

  arma::mat EE(T, N, arma::fill::zeros) ;
  arma::mat FE(T, N, arma::fill::zeros) ; // store overall effects

  // arma::mat F1(1, r, arma::fill::zeros) ;
  // arma::mat F2(T - 1, r, arma::fill::zeros) ;
  // arma::mat L1(1, r, arma::fill::zeros) ;
  // arma::mat L2(N - 1, r, arma::fill::zeros) ;

  arma::mat F1 ;
  arma::mat F2 ;
  arma::mat L1 ;
  arma::mat L2 ;

  arma::mat EE1 ;
  arma::mat EE2 ;

  // step1. update lambda 
  // remove additive fixed effects
  if (r > 0 || force == 1 || force == 3) {
    EE = E ;
    if (force == 2 || force == 3) {
      EE = EE - repmat(xi_old, 1, N) ;
    }
    if (r > 0) {
      if (force == 0) {
        // no additional restrictions
        // L = arma::pinv(F_old) * EE ;
        // L = L.t() ;
        L = EE.t() * (arma::pinv(F_old)).t()  ;
      } else {
        EE1 = EE.col(0) ;
        EE2 = EE.cols(1, N - 1) ;
        if (force == 1) {
          // alpha_1 = 0, f_1 = 0
          // L1 = arma::pinv(F_old) * EE1 ;
          // L1 = L1.t() ;
          L1 = EE1.t() * (arma::pinv(F_old)).t() ;
          // intercept is the unit fixed effect
          F_old.insert_cols(0, et) ;
          // L2 = arma::pinv(F_old) * EE2 ;
          // L2 = L2.t() ;
          L2 = EE2.t() * (arma::pinv(F_old)).t() ;
          // unit fixed effect 
          alpha = L2.col(0) ;
          // Rcpp::Rcout << "OK1" << std::endl;
          alpha.insert_rows(0, u0) ;
          // remove unit fe from loadings
          L2.shed_col(0) ;
          // Rcpp::Rcout << "OK2" << std::endl;
          L2.insert_rows(0, L1) ;
          L = L2 ;
        }
        else if (force == 2) {
          // xi_1 = 0, lambda_1 = 0
          // L2 = arma::pinv(F_old) * EE2 ;
          // L2 = L2.t() ;
          L2 = EE2.t() * (arma::pinv(F_old)).t() ;
          L2.insert_rows(0, f0) ;
          L = L2 ;
        }
        else if (force == 3) {
          // alpha_1 = 0, xi_1 = 0, lambda_1 = 0, f_1 = 0
          // intercept is the unit fixed effect
          F_old.insert_cols(0, et) ;
          // L2 = arma::pinv(F_old) * EE2 ;
          // L2 = L2.t() ;
          L2 = EE2.t() * (arma::pinv(F_old)).t() ;
          // unit fixed effect 
          alpha = L2.col(0) ;
          alpha.insert_rows(0, u0) ;
          // remove unit fe from loadings
          L2.shed_col(0) ;
          L2.insert_rows(0, f0) ;
          L = L2 ;  
        }
      }
    } else {
      if (force == 1 || force == 3) {
        // alpha_1 = 0
        alpha  =  mean(EE, 0).t() ;
        alpha(0, 0) = 0 ;
      }
    }
  }

  // step2. update factor 
  // remove updated alpha
  if (r > 0 || force == 1 || force == 3) {
    EE = E ;
    if (force == 1 || force == 3) {
      EE = EE - repmat(alpha.t(), T, 1) ;
    } 
    if (r > 0) {
      if (force == 0) {
        // F = arma::pinv(L) * EE.t() ;
        // F = F.t() ;
        F = EE * (arma::pinv(L)).t() ;
      } else {
        EE1 = EE.row(0) ;
        EE2 = EE.rows(1, T - 1) ;
        if (force == 1) {
          // alpha_1 = 0, f_1 = 0
          // F2 = arma::pinv(L) * EE2.t() ;
          // F2 = F2.t() ;
          F2 = EE2 * (arma::pinv(L)).t() ;
          // Rcpp::Rcout << "OK3" << std::endl;
          F2.insert_rows(0, f0) ;
          F = F2 ;
        }
        else if (force == 2) {
          // xi_1 = 0
          // F1 = arma::pinv(L) * EE1.t() ;
          // F1 = F1.t() ;
          F1 = EE1 * (arma::pinv(L)).t() ;
          L.insert_cols(0, en) ;
          // F2 = arma::pinv(L) * EE2.t() ;
          // F2 = F2.t() ;
          F2 = EE2 * (arma::pinv(L)).t() ;
          // time fixed effect 
          xi = F2.col(0) ;
          xi.insert_rows(0, u0) ;
          // remove time fe from factor
          F2.shed_col(0) ;
          F2.insert_rows(0, F1) ;
          F = F2 ;
          L.shed_col(0) ;
        }
        else if (force == 3) {
          // alpha_1 = 0, xi_1 = 0, lambda_1 = 0, f_1 = 0
          // intercept is the time fixed effect
          L.insert_cols(0, en) ;
          // F2 = arma::pinv(L) * EE2.t() ;
          // F2 = F2.t() ;
          F2 = EE2 * (arma::pinv(L)).t() ;
          // time fixed effect 
          xi = F2.col(0) ;
          xi.insert_rows(0, u0) ;
          // remove time fe from factor
          F2.shed_col(0) ;
          F2.insert_rows(0, f0) ;
          F = F2 ;
          L.shed_col(0) ;  
        }
      }
    } else {
      if (force == 2 || force == 3) {
        // xi_1 = 0
        xi  =  mean(EE, 1) ;
        xi(0, 0) = 0 ;
      }
    }
  }

  List result ;
  if (force == 1 || force == 3) {
    result["alpha"] = alpha ;
    FE = FE + repmat(alpha.t(), T, 1) ;
  }
  if (force == 2 || force == 3) {
    result["xi"] = xi ;
    FE = FE + repmat(xi, 1, N) ;
  }
  if (r > 0) {
    List qrf ;
    qrf = qr_factor(F, L) ;
    result["factor"] = as<arma::mat>(qrf["factor"]) ;
    result["lambda"] = as<arma::mat>(qrf["lambda"]) ;
    FE = FE + as<arma::mat>(qrf["FE"]) ;
  }
  result["FE"] = FE ;
  return(result) ;

}

/* factor analysis: mu add ife, unbalanced */
// [[Rcpp::export]]
List fe_ub (const arma::mat& E,
            const arma::mat& I,
            const arma::mat& F_old,
            const arma::mat& xi_old,
            int force,
            int r
           ) {

  int T = E.n_rows ; 
  int N = E.n_cols ;

  arma::mat F(T, r, arma::fill::zeros) ;
  arma::mat L(N, r, arma::fill::zeros) ;
  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;

  arma::mat EE(T, N, arma::fill::zeros) ;
  arma::mat FE(T, N, arma::fill::zeros) ; // store overall effects

  // update lambda 
  if (r > 0 || force == 1 || force == 3) {
    EE = E ;
    // s1 remove additive fe
    if (force == 2 || force == 3) {
      EE = EE - repmat(xi_old, 1, N) ;
    }
    EE = FE_adj(EE, I) ;

    List upl ;
    upl = l_ub (EE, F_old, I, r, force) ;
    if (r > 0) {
      L = as<arma::mat>(upl["L"]) ;
    }
    if (force == 1 || force == 3) {
      alpha = as<arma::mat>(upl["alpha"]) ;
    }
  }

  // update factor  
  if (r > 0 || force == 2 || force == 3) {
    EE = E ;
    // s1 remove updated additive unit fe
    if (force == 1 || force == 3) {
      // alpha = as<arma::mat>(upl["alpha"]) ;
      EE = EE - repmat(alpha.t(), T, 1) ;
    } 
    EE = FE_adj(EE, I) ;

    List upf ;
    upf = f_ub(EE, L, I, r, force) ;
    if (r > 0) {
      F = as<arma::mat>(upf["F"]) ;
    }
    if (force == 2 || force == 3) {
      xi = as<arma::mat>(upf["xi"]) ;
    }
  }

  List result ;

  if (force == 1 || force == 3) {
    FE = FE + repmat(alpha.t(), T, 1) ;
    result["alpha"] = alpha ;
    // FE = FE + repmat(alpha.t(), T, 1) ;
  }
  if (force == 2 || force == 3) {
    FE = FE + repmat(xi, 1, N) ;
    result["xi"] = xi ;
    // FE = FE + repmat(xi, 1, N) ;
  }
  if (r > 0) {
    List qrf ;
    qrf = qr_factor(F, L) ;
    result["factor"] = as<arma::mat>(qrf["factor"]) ;
    result["lambda"] = as<arma::mat>(qrf["lambda"]) ;
    FE = FE + as<arma::mat>(qrf["FE"]) ;
  }
  result["FE"] = FE ;
  return(result) ;
}
