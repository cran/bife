#define ARMA_NO_DEBUG


#include <RcppArmadillo.h>
#include "link.h"
#include "demeaning.h"
#include "biascorr_apeff.h"
#include "offset.h"


// Uncorrected average partial effects
// [[Rcpp::export(name = ".apeff")]]
arma::colvec apeff(const arma::colvec &discrete,
                   const arma::colvec &beta,
                   const arma::colvec &alpha,
                   const arma::mat &X,
                   const arma::uvec &index_id,
                   const unsigned int model) {
 
 
 const arma::uword D = X.n_cols;
 const arma::uword NT = X.n_rows;
 
 arma::mat ipeff(NT, D);
 for (arma::uword i = 0 ; i < D ; ++i) {
  
  if (discrete(i) == 1) {
   
   arma::mat X_one = X;
   arma::mat X_zero = X;
   X_one.col(i).ones();
   X_zero.col(i).zeros();
   
   switch (model) {
   
   case 0 :
    ipeff.col(i) = logistic_cdf(X_one * beta + alpha(index_id)) - logistic_cdf(X_zero * beta + alpha(index_id));
    break;
   case 1 :
    ipeff.col(i) = normal_cdf(X_one * beta + alpha(index_id)) - normal_cdf(X_zero * beta + alpha(index_id));
    break;
   }
  } else {
   
   switch (model) {
   
   case 0 :
    ipeff.col(i) = beta(i) * logistic_pdf(X * beta + alpha(index_id));
    break;
   case 1 :
    ipeff.col(i) = beta(i) * normal_pdf(X * beta + alpha(index_id));
    break;
   }
  }
 }
 
 
 return arma::mean(ipeff).t();
}


// Corrected average partial effects
// [[Rcpp::export(name = ".apeff.corr")]]
arma::mat apeff_tilde(const arma::colvec &discrete,
                      const arma::colvec &beta,
                      const arma::colvec &alpha,
                      const arma::colvec &beta_start,
                      const arma::colvec &beta_tilde,
                      const arma::colvec &alpha_tilde,
                      const arma::colvec &y,
                      const arma::mat &X,
                      const arma::colvec &t,
                      const arma::uvec &index_id,
                      const arma::colvec &T_vector,
                      const unsigned int model,
                      const unsigned int bias_corr,   // 0 = semi corrected, 1 = analytically, 2 = jackknife
                      const unsigned int iter_max1,  
                      const double tolerance1,
                      const unsigned int iter_max2,  
                      const double tolerance2) {    
 
 
 const arma::uword D = beta.n_rows;
 
 
 /* --------------------------------------------------------------------------------------------------------
 *  Bias Correction 
 *
 *  - No bias correction / parial effects semi corrected   - bias_corr = 0
 *  - Analytical bias correction - Hahn and Newey (2004)   - bias_corr = 1
 *  - Jackknife bias correction - Hahn and Newey (2004)    - bias_corr = 2 
 * --------------------------------------------------------------------------------------------------------
 */ 
 
 arma::mat avgpeff(D, 2);
 avgpeff.col(0) = apeff(discrete, beta, alpha, X, index_id, model);
 switch (bias_corr) {
 
 case 0 :
  avgpeff.col(1) = apeff(discrete, beta_tilde, alpha_tilde, X, index_id, model);
  break;
 case 1 :
  avgpeff.col(1) = apeff(discrete, beta_tilde, alpha_tilde, X, index_id, model) - biascorr_apeff(discrete, beta_tilde, alpha_tilde, y, X, index_id, T_vector, model);
  break;
 case 2 :
  const arma::uword N = T_vector.n_rows;
  arma::colvec y_mean(N);
  arma::mat X_mean(N, D);
  arma::uword start = 0;
  for (arma::uword i = 0 ; i < N ; ++i) {
   
   const arma::uword end = T_vector(i) + start - 1;
   y_mean(i) = arma::mean(y.subvec(start, end));
   X_mean.row(i) = arma::mean(X.rows(start, end));
   start = end + 1;
  }
  const arma::uword T = t.max();
  const arma::colvec T_jack = T_vector - 1;
  arma::colvec se_beta_jack(D);
  arma::colvec se_alpha_jack(N);
  arma::mat H_jack(D, D);
  unsigned int iter_jack;
  bool conv_jack;
  
  arma::mat avgpeff_jack(D, T);
  for (arma::uword i = 0 ; i < T ; ++i) {
   
   Rcpp::checkUserInterrupt();
   
   const arma::uvec index_jack = arma::find(t != i + 1);
   const arma::colvec y_jack = y(index_jack);
   const arma::mat X_jack = X.rows(index_jack);
   const arma::uvec index_id_jack = index_id(index_jack);
   
   arma::colvec beta_jack = beta_start;
   arma::colvec alpha_jack(N);
   switch (model) {
   
   case 0 : 
    alpha_jack = logistic_inv(y_mean) - X_mean * beta_jack;
    break;
   case 1 :
    alpha_jack = normal_inv(y_mean) - X_mean * beta_jack;
    break;
   }
   
   demeaning(beta_jack, alpha_jack, se_beta_jack, se_alpha_jack, H_jack, iter_jack, conv_jack, y_jack, X_jack, index_id_jack, T_jack, model, iter_max1, tolerance1);
   alpha_jack = glm_offset(iter_jack, conv_jack, beta_jack, y_jack, X_jack, index_id_jack, T_jack, model, iter_max2, tolerance2);
   
   avgpeff_jack.col(i) = apeff(discrete, beta_jack, alpha_jack, X_jack, index_id_jack, model);
  }
  const double mean_T = arma::mean(T_vector);
  avgpeff.col(1) = mean_T * apeff(discrete, beta, alpha, X, index_id, model) - (mean_T - 1) * arma::sum(avgpeff_jack, 1) / mean_T;
  break;
 }
 
 
 return avgpeff;
}
