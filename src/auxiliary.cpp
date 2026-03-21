#include "fect.h"

// auxiliary functions


/* cross product */
arma::mat crossprod (const arma::mat& x, const arma::mat& y) {
  return(x.t() * y);
}

/* Expectation :E if Iij==0, Eij=FEij */
arma::mat E_adj (const arma::mat& E, const arma::mat& FE,
                 const arma::mat& I) {
  int T = E.n_rows ;
  int N = E.n_cols ;
  arma::mat EE = E ;
  for (int i = 0; i < T; i++) {
    for (int j = 0; j < N; j++) {
      if (I(i, j) == 0) {
        EE(i, j) = FE(i, j) ;
      }
    }
  }
  return(EE) ;
}

/* Expectation Step :E if Iij==0, Eij=FEij with Weights */
// [[Rcpp::export]]
arma::mat wE_adj (const arma::mat& E, const arma::mat& FE, const arma::mat& W,
                  const arma::mat& I) {
  int T = E.n_rows ;
  int N = E.n_cols ;
  arma::mat EE = E ;
  for (int i = 0; i < T; i++) {
    for (int j = 0; j < N; j++) {
      if (I(i, j) == 0) {
        EE(i, j) = FE(i, j) ;
      }
      else{
        EE(i, j) = (1-W(i,j))*FE(i, j) + W(i,j)*E(i, j) ;
      }
    }
  }
  return(EE) ;
}


/* reset FEij=0 if Iij==0 */
arma::mat FE_adj (const arma::mat& FE, const arma::mat& I) {
  int T = FE.n_rows ;
  int N = FE.n_cols ;
  arma::mat FEE = FE ;
  for (int i = 0; i < T; i++) {
    for (int j = 0; j < N; j++) {
      if (I(i, j) == 0) {
        FEE(i, j) = 0 ;
      }
    }
  }
  return(FEE) ;
}

/* drop values if Iij == 1 */
arma::mat FE_missing (const arma::mat& FE, const arma::mat& I) {
  int T = FE.n_rows ;
  int N = FE.n_cols ;
  arma::mat FEE = FE ;
  for (int i = 0; i < T; i++) {
    for (int j = 0; j < N; j++) {
      if (I(i, j) == 1) {
        FEE(i, j) = 0 ;
      }
    }
  }
  return(FEE) ;
}
        /* ------------------------------------------ */
/* ------------------- 2. Only for Probit ----------------------- */
        /* ------------------------------------------ */

arma::mat subr (const arma::mat& X, const arma::mat& ind) {
  int N = X.n_cols ;
  int t = ind.n_rows ;
  arma::mat subX(t, N, arma::fill::zeros) ;

  for (int i = 0; i < t; i++) {
    subX.row(i) = X.row(ind(i,0)) ;
  }
  return(subX) ;
}

/* compute M in the EM algorithm for probit model */
arma::mat M_gen (const arma::mat& Y_fit, const arma::mat& Y) {
  int T = Y.n_rows;
  int N = Y.n_cols;
  arma::mat M(T, N, arma::fill::zeros);
  for(int i=0; i<T; i++) {
    for(int j=0; j<N; j++) {
      if(Y(i,j)==1){
        M(i,j)=R::dnorm(-Y_fit(i,j),0.0,1.0,false)
        /(1-R::pnorm(-Y_fit(i,j),0.0,1.0,true,false));
      }
      else if(Y(i,j)==0){
        M(i,j)=-R::dnorm(-Y_fit(i,j),0.0,1.0,false)
        /R::pnorm(-Y_fit(i,j),0.0,1.0,true,false);
      }
    }
  }
  return(M);
}

/* compute M in the EM algorithm for probit model for ub data */
arma::mat M_gen_ub (const arma::mat& Y_fit, const arma::mat& Y, const arma::mat& I) {
  int T = Y.n_rows;
  int N = Y.n_cols;
  arma::mat M(T, N, arma::fill::zeros);
  for(int i=0; i<T; i++) {
    for(int j=0; j<N; j++) {
      if(I(i,j)!=0){
        if(Y(i,j)==1){
          M(i,j)=R::dnorm(-Y_fit(i,j),0.0,1.0,false)
          /(1-R::pnorm(-Y_fit(i,j),0.0,1.0,true,false));
        }
        else if(Y(i,j)==0){
          M(i,j)=-R::dnorm(-Y_fit(i,j),0.0,1.0,false)
          /R::pnorm(-Y_fit(i,j),0.0,1.0,true,false);
        }
      }  
    }
  }
  return(M);
}

/* calculate alpha in px */
double alpha_hat (const arma::mat& res, const arma::mat& v) {
  int T = res.n_rows;
  int N = res.n_cols;
  arma::mat vhat(T, N, arma::fill::zeros);
  for (int i = 0; i < T; i++) {
    for (int j = 0; j < N; j++) {
      vhat(i,j) = pow(res(i,j),2) + v(i,j);
    }
  }
  double alhat2 = arma::accu(vhat)/(N*T);
  double alhat = sqrt(alhat2);
  return(alhat);
}

/* calculate S in mopx */
double S (double a, double b, double w) {
  double out;
  if(b>=a){
    out = (1+w)*b-w*a;
  }
  else{
    out = pow(b,(1+w))*pow(a,(-w));
  }
  return(out);
}

/* calculate gamma_hat in mopx */
double gamma_hat (const arma::mat& res, const arma::mat& V) {
  int T = res.n_rows;
  int N = res.n_cols;
  arma::mat gam(T,N, arma::fill::zeros);
  double out;
  for(int i=0; i<T; i++){
    for(int j=0; j<N; j++){
      gam(i,j)=pow(res(i,j),2)+V(i,j);
    }
  }
  out = (N*T)/arma::accu(gam);
  return(out);
}

/* calculate gamma_hat in mopx; unbalanced */
double gamma_hat_ub (const arma::mat& res, const arma::mat& V, const arma::mat& I) {
  int T = res.n_rows;
  int N = res.n_cols;
  arma::mat gam(T, N, arma::fill::zeros);
  double out;
  for(int i=0; i<T; i++){
    for(int j=0; j<N; j++){
      if(I(i,j)==1){
        gam(i,j)=pow(res(i,j),2)+V(i,j);
      }
    }
  }
  out = arma::accu(I)/arma::accu(gam);
  return(out);
}


/* calculate v in px */
arma::mat V (const arma::mat& Y, const arma::mat& Y_fit) {
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  arma::mat Vij(T, N, arma::fill::zeros) ;
  for(int i=0; i<T; i++){
    for(int j=0; j<N; j++){
      if(Y(i,j)==1){
        Vij(i,j)=1-(Y_fit(i,j)+R::dnorm(Y_fit(i,j),0.0,1.0,false)
          /R::pnorm(Y_fit(i,j),0.0,1.0,true,false))*
            (R::dnorm(Y_fit(i,j),0.0,1.0,false)/
             R::pnorm(Y_fit(i,j),0.0,1.0,true,false)) ;
      }
      else if(Y(i,j)==0){
        Vij(i,j)=1-(-Y_fit(i,j)+R::dnorm(Y_fit(i,j),0.0,1.0,false)
          /R::pnorm(-Y_fit(i,j),0.0,1.0,true,false))*
            (R::dnorm(Y_fit(i,j),0.0,1.0,false)/
             R::pnorm(-Y_fit(i,j),0.0,1.0,true,false)) ;
      }
    }
  }
  return(Vij) ;
}

/* calculate v in px; unbalanced */
arma::mat V_ub (const arma::mat& Y, const arma::mat& Y_fit, const arma::mat& I) {
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  arma::mat Vij(T, N, arma::fill::ones) ;
  for(int i=0; i<T; i++){
    for(int j=0; j<N; j++){
      if(I(i,j)==1){
        if(Y(i,j)==1){
          Vij(i,j)=1-(Y_fit(i,j)+R::dnorm(Y_fit(i,j),0.0,1.0,false)
              /R::pnorm(Y_fit(i,j),0.0,1.0,true,false))*
              (R::dnorm(Y_fit(i,j),0.0,1.0,false)/
              R::pnorm(Y_fit(i,j),0.0,1.0,true,false)) ;
        }
        else if(Y(i,j)==0){
          Vij(i,j)=1-(-Y_fit(i,j)+R::dnorm(Y_fit(i,j),0.0,1.0,false)
            /R::pnorm(-Y_fit(i,j),0.0,1.0,true,false))*
              (R::dnorm(Y_fit(i,j),0.0,1.0,false)/
               R::pnorm(-Y_fit(i,j),0.0,1.0,true,false)) ;
        }
      }
    }
  }
  return(Vij) ;
}
        /* ------------------------------------------ */
/* ------------------ 3. Functions exported  ----------------------- */
        /* ------------------------------------------ */

/* probit: log likelihood: sum */
// [[Rcpp::export]]
double loglh (const arma::mat& Y_fit, const arma::mat& Y) {
  int T = Y.n_rows;
  int N = Y.n_cols;
  arma::mat lh(T, N, arma::fill::zeros);
  for(int i=0; i<T; i++){
    for(int j=0; j<N; j++){
      if(Y(i,j)==1){
        lh(i,j)=log(R::pnorm(Y_fit(i,j),0.0,1.0,true,false));
      }
      else if(Y(i,j)==0){
          lh(i,j)=log(1-R::pnorm(Y_fit(i,j),0.0,1.0,true,false));
      }
    }
  }
  double llh = arma::accu(lh);
  return(llh); 
}

/* probit: log likelihood: sum unbalanced case */
// [[Rcpp::export]]
double loglh_ub (const arma::mat& Y_fit, const arma::mat& Y, const arma::mat& I) {
  int T = Y.n_rows;
  int N = Y.n_cols;
  arma::mat lh(T, N, arma::fill::zeros);
  for(int i=0; i<T; i++){
    for(int j=0; j<N; j++){
      if(I(i,j)==1){
        if(Y(i,j)==1){
          lh(i,j)=log(R::pnorm(Y_fit(i,j),0.0,1.0,true,false));
        }
        else if(Y(i,j)==0){
          lh(i,j)=log(1-R::pnorm(Y_fit(i,j),0.0,1.0,true,false));
        }
      }
    }
  }
  double llh = arma::accu(lh);
  return(llh); 
} 

/* adjust unbalanced data */
// [[Rcpp::export]]
arma::mat data_ub_adj (const arma::mat& I_data, const arma::mat& data) {
  int count = I_data.n_rows ;
  //int total = data.n_rows ;
  int nov = data.n_cols ;
  arma::mat data_adj(count,nov) ;
  data_adj.fill(arma::datum::nan) ;
  int j = 0;
  for(int i=0; i<count; i++){
    if(I_data(i,0)==1){
      data_adj.row(i) = data.row(j);
      j++;
    }
  }
  return(data_adj);
}

/* Three dimensional matrix inverse */
// [[Rcpp::export]]
arma::mat XXinv (const arma::cube& X) { 
  int p = X.n_slices ;
  arma::mat xx(p, p) ;
  for (int k = 0; k < p; k++) {
    for (int m = 0; m < p; m++) {
      if (k > m) {
        xx(k, m) = xx(m, k);
      }
      else {
        xx(k, m) = trace(crossprod(X.slice(k), X.slice(m))) ;
      }
    }
  } 
  return(inv(xx)) ;
}

/* weighted inverse*/
/* X: T*N*p; W:T*N */
// [[Rcpp::export]]
arma::mat wXXinv (const arma::cube& X, const arma::mat& w) { 
  int p = X.n_slices ;
  arma::mat w_sr = pow(w,0.5) ;
  arma::mat xx(p, p) ;
  for (int k = 0; k < p; k++) {
    for (int m = 0; m < p; m++) {
      if (k > m) {
        xx(k, m) = xx(m, k);
      }
      else {
        xx(k, m) = trace(crossprod(w_sr%X.slice(k), w_sr%X.slice(m))) ;
      }
    }
  }  
  return(inv_sympd(xx)) ;
}


/* Obtain beta given interactive fe */
// [[Rcpp::export]]
arma::mat panel_beta (const arma::cube& X, const arma::mat& xxinv,
                      const arma::mat& Y, const arma::mat& FE) {
  int p = X.n_slices ; 
  arma::mat xy(p, 1, arma::fill::zeros) ;
  for (int k = 0; k < p; k++) {
    xy(k) = trace(crossprod(X.slice(k), (Y - FE))) ;
  }
  return(xxinv * xy);
}

/* Obtain beta given interactive fe using weighted regression */
// [[Rcpp::export]]
arma::mat wpanel_beta (const arma::cube& X, const arma::mat& xwxinv, const arma::mat& w,
                       const arma::mat& Y, const arma::mat& FE) {
  int p = X.n_slices ; 
  arma::mat xwy(p, 1, arma::fill::zeros) ;
  arma::mat w_sr = pow(w,0.5) ;
  for (int k = 0; k < p; k++) {
    xwy(k) = trace(crossprod(w_sr%X.slice(k), w_sr%(Y - FE))) ;
  }
  return(xwxinv * xwy);
}


/* Obtain OLS panel estimate */
// [[Rcpp::export]]
arma::mat panel_est (const arma::cube& X, const arma::mat& Y, const arma::mat& MF) {
  int p = X.n_slices ;
  arma::mat xx(p, p, arma::fill::zeros);
  arma::mat xy(p, 1, arma::fill::zeros);
  if (p==1) {
    xx(0,0) = trace(X.slice(0).t() * MF * X.slice(0)) ;
    xy(0,0) = trace(X.slice(0).t() * MF * Y) ;
  }
  if (p>1) {
    for (int k = 0; k < p; k++) {
      arma::mat MFX1 = MF * X.slice(k) ;
      xy(k, 0) = trace(crossprod(MFX1, Y)) ;
      for (int m = k; m < p; m++) {
        arma::mat MFX2 = MF * X.slice(m) ;
        xx(k, m) = trace(crossprod(MFX1, MFX2)) ;
        if (k < m) {
          xx(m, k) = xx(k, m) ;
        }
      }
    }
  }
  return(xx.i() * xy) ;
}
