#ifndef __FECT_H__
#define __FECT_H__


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp ;

arma::mat crossprod (const arma::mat& x, const arma::mat& y);
arma::mat E_adj (const arma::mat& E, const arma::mat& FE, const arma::mat& I);
arma::mat wE_adj (const arma::mat& E, const arma::mat& FE, const arma::mat& W, const arma::mat& I);
arma::mat FE_adj (const arma::mat& FE, const arma::mat& I);
arma::mat FE_missing (const arma::mat& FE, const arma::mat& I);
arma::mat subr (const arma::mat& X, const arma::mat& ind);
arma::mat M_gen (const arma::mat& Y_fit, const arma::mat& Y);
arma::mat M_gen_ub (const arma::mat& Y_fit, const arma::mat& Y, const arma::mat& I);
double alpha_hat (const arma::mat& res, const arma::mat& v);
double S (double a, double b, double w);
double gamma_hat (const arma::mat& res, const arma::mat& V);
double gamma_hat_ub (const arma::mat& res, const arma::mat& V, const arma::mat& I);
arma::mat V (const arma::mat& Y, const arma::mat& Y_fit);
arma::mat V_ub (const arma::mat& Y, const arma::mat& Y_fit, const arma::mat& I);
double loglh(const arma::mat& Y_fit, const arma::mat& Y);
double loglh_ub(const arma::mat& Y_fit, const arma::mat& Y, const arma::mat& I);
arma::mat data_ub_adj(const arma::mat& I_data, const arma::mat& data);
arma::mat XXinv(const arma::cube& X);
arma::mat wXXinv(const arma::cube& X, const arma::mat& w); 
arma::mat panel_beta(const arma::cube& X, const arma::mat& xxinv, const arma::mat& Y, const arma::mat& FE);
arma::mat wpanel_beta (const arma::cube& X, const arma::mat& xwxinv, const arma::mat& w, const arma::mat& Y, const arma::mat& FE);
List Y_demean(const arma::mat& Y, int force);
List Y_wdemean(const arma::mat& Y, const arma::mat& W, int force);
List fe_add(const arma::mat& alpha_Y, const arma::mat& xi_Y, double mu_Y, int T, int N, int force);
arma::mat panel_est(const arma::cube& X, const arma::mat& Y, const arma::mat& MF);
List panel_factor(const arma::mat& E, int r);
arma::mat panel_FE(const arma::mat& E, double lambda, int hard);
List ife(const arma::mat& E, int force, int mc, int r, int hard, double lambda);
List fe_ad_iter(const arma::mat& Y, const arma::mat& Y0, const arma::mat& I, const arma::mat& W, int force, double tolerate, int max_iter);
List fe_ad_covar_iter(const arma::cube& XX, const arma::mat& xxinv, const arma::mat& Y, const arma::mat& Y0, const arma::mat& I, const arma::mat& beta0, const arma::mat& W, int force, double tolerate, int max_iter);
List fe_ad_inter_iter(const arma::mat& Y, const arma::mat& Y0, const arma::mat& I, const arma::mat& W, int force, int mc, int r, int hard, double lambda, double tolerate, int max_iter);
List fe_ad_inter_covar_iter(const arma::cube& XX, const arma::mat& xxinv, const arma::mat& Y, const arma::mat& Y0, const arma::mat& I, const arma::mat& W, const arma::mat& beta0, int force, int mc, int r, int hard, double lambda, double tolerate, int max_iter);
List beta_iter(const arma::cube& X, const arma::mat& xxinv, const arma::mat& Y, int r, double tolerate, const arma::mat& beta0, int max_iter);
List qr_factor(const arma::mat& F, const arma::mat& L);
arma::mat IND(const arma::mat& I);
List subfe(const arma::mat& Y, const arma::mat& X, const arma::mat& I, int intercept);
List l_ub(const arma::mat& Y, const arma::mat& F, const arma::mat& I, int r, int force);
List f_ub(const arma::mat& Y, const arma::mat& L, const arma::mat& I, int r, int force);
List fe(const arma::mat& E, const arma::mat& F_old, const arma::mat& xi_old, int force, int r);
List fe_ub(const arma::mat& E, const arma::mat& I, const arma::mat& F_old, const arma::mat& xi_old, int force, int r);
List inter_fe(const arma::mat& Y, const arma::cube& X, int r, int force, const arma::mat& beta0, double tol, int max_iter);
List inter_fe_ub(const arma::mat& Y, const arma::mat& Y0, const arma::cube& X, const arma::mat& I, const arma::mat& W, const arma::mat& beta0, int r, int force, double tol, int max_iter);
List inter_fe_mc(const arma::mat& Y, const arma::mat& Y0, const arma::cube& X, const arma::mat& I, const arma::mat& W, const arma::mat& beta0, int r, double lambda, int force, double tol, int max_iter);
List inter_fe_d(const arma::mat& Y, const arma::mat& Y_fit0, const arma::mat& FE0, const arma::cube& X, int r, int force, int mniter, double w, double tol);
List inter_fe_d_ub(const arma::mat& Y, const arma::mat& Y_fit0, const arma::mat& FE0, const arma::cube& X, const arma::mat& I, int r, int force, int mniter, double w, double tol);
List inter_fe_d_qr(const arma::mat& Y, const arma::mat& Y_fit0, const arma::mat& FE0, const arma::mat& factor0, const arma::mat& xi0, const arma::cube& X, int r, int force, int mniter, double w, double tol);
List inter_fe_d_qr_ub(const arma::mat& Y, const arma::mat& Y_fit0, const arma::mat& FE0, const arma::mat& factor0, const arma::mat& xi0, const arma::cube& X, const arma::mat& I, int r, int force, int mniter, double w, double tol);
List cfe_iter(const arma::cube& XX, const arma::mat& xxinv, const arma::cube& X_extra_FE, const arma::cube& X_Z, const arma::cube& X_Q, const arma::cube& X_gamma, const arma::cube& X_kappa, Rcpp::List Zgamma_id, Rcpp::List kappaQ_id, const arma::mat& Y, const arma::mat& Y0, const arma::mat& I, const arma::mat& W, const arma::mat& beta0, int force, int r, double tolerate, int max_iter);
List complex_fe_ub(const arma::mat& Y, const arma::mat& Y0, const arma::cube& X_covariates, const arma::cube& X_extra_FE, const arma::cube& X_Z, const arma::cube& X_Q, const arma::cube& X_gamma, const arma::cube& X_kappa, Rcpp::List Zgamma_id, Rcpp::List kappaQ_id, const arma::mat& I, const arma::mat& W, const arma::mat& beta0, int r, int force, double tol, int max_iter);
#endif
