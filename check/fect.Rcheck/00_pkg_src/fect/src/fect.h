#ifndef __FECT_H__
#define __FECT_H__


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
using namespace Rcpp ;

arma::mat crossprod (arma::mat x, arma::mat y);
arma::mat E_adj (arma::mat E, arma::mat FE, arma::mat I);
arma::mat wE_adj (arma::mat E, arma::mat FE, arma::mat W, arma::mat I);
arma::mat FE_adj (arma::mat FE, arma::mat I);
arma::mat FE_missing (arma::mat FE, arma::mat I);
arma::mat subr (arma::mat X, arma::mat ind);
arma::mat M_gen (arma::mat Y_fit, arma::mat Y);
arma::mat M_gen_ub (arma::mat Y_fit, arma::mat Y, arma::mat I);
double alpha_hat (arma::mat res, arma::mat v);
double S (double a, double b, double w);
double gamma_hat (arma::mat res, arma::mat V);
double gamma_hat_ub (arma::mat res, arma::mat V, arma::mat I);
arma::mat V (arma::mat Y, arma::mat Y_fit);
arma::mat V_ub (arma::mat Y, arma::mat Y_fit, arma::mat I);
double loglh(arma::mat Y_fit, arma::mat Y);
double loglh_ub(arma::mat Y_fit, arma::mat Y, arma::mat I);
arma::mat data_ub_adj(arma::mat I_data, arma::mat data);
arma::mat XXinv(arma::cube X);
arma::mat wXXinv(arma::cube X, arma::mat w); 
arma::mat panel_beta(arma::cube X, arma::mat xxinv, arma::mat Y, arma::mat FE);
arma::mat wpanel_beta (arma::cube X, arma::mat xwxinv, arma::mat w, arma::mat Y, arma::mat FE);
List Y_demean(arma::mat Y, int force);
List Y_wdemean(arma::mat Y, arma::mat W, int force);
List fe_add(arma::mat alpha_Y, arma::mat xi_Y, double mu_Y, int T, int N, int force);
arma::mat panel_est(arma::cube X, arma::mat Y, arma::mat MF);
List panel_factor(arma::mat E, int r);
arma::mat panel_FE(arma::mat E, double lambda, int hard);
List ife(arma::mat E, int force, int mc, int r, int hard, double lambda);
List fe_ad_iter(arma::mat Y, arma::mat Y0, arma::mat I, arma::mat W, int force, double tolerate, int max_iter);
List fe_ad_covar_iter(arma::cube XX, arma::mat xxinv, arma::mat Y, arma::mat Y0, arma::mat I, arma::mat beta0, arma::mat W, int force, double tolerate, int max_iter);
List fe_ad_inter_iter(arma::mat Y, arma::mat Y0, arma::mat I, arma::mat W, int force, int mc, int r, int hard, double lambda, double tolerate, int max_iter);
List fe_ad_inter_covar_iter(arma::cube XX, arma::mat xxinv, arma::mat Y, arma::mat Y0, arma::mat I, arma::mat W, arma::mat beta0, int force, int mc, int r, int hard, double lambda, double tolerate, int max_iter);
List beta_iter(arma::cube X, arma::mat xxinv, arma::mat Y, int r, double tolerate, arma::mat beta0, int max_iter);
List qr_factor(arma::mat F, arma::mat L);
arma::mat IND(arma::mat I);
List subfe(arma::mat Y, arma::mat X, arma::mat I, int intercept);
List l_ub(arma::mat Y, arma::mat F, arma::mat I, int r, int force);
List f_ub(arma::mat Y, arma::mat L, arma::mat I, int r, int force);
List fe(arma::mat E, arma::mat F_old, arma::mat xi_old, int force, int r);
List fe_ub(arma::mat E, arma::mat I, arma::mat F_old, arma::mat xi_old, int force, int r);
List inter_fe(arma::mat Y, arma::cube X, int r, int force, arma::mat beta0, double tol, int max_iter);
List inter_fe_ub(arma::mat Y, arma::mat Y0, arma::cube X, arma::mat I, arma::mat W, arma::mat beta0, int r, int force, double tol, int max_iter);
List inter_fe_mc(arma::mat Y, arma::mat Y0, arma::cube X, arma::mat I, arma::mat W, arma::mat beta0, int r, double lambda, int force, double tol, int max_iter);
List inter_fe_d(arma::mat Y, arma::mat Y_fit0, arma::mat FE0, arma::cube X, int r, int force, int mniter, double w, double tol);
List inter_fe_d_ub(arma::mat Y, arma::mat Y_fit0, arma::mat FE0, arma::cube X, arma::mat I, int r, int force, int mniter, double w, double tol);
List inter_fe_d_qr(arma::mat Y, arma::mat Y_fit0, arma::mat FE0, arma::mat factor0, arma::mat xi0, arma::cube X, int r, int force, int mniter, double w, double tol);
List inter_fe_d_qr_ub(arma::mat Y, arma::mat Y_fit0, arma::mat FE0, arma::mat factor0, arma::mat xi0, arma::cube X, arma::mat I, int r, int force, int mniter, double w, double tol);

#endif
