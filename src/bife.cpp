#define ARMA_NO_DEBUG


#include <RcppArmadillo.h>
#include "link.h"
#include "demeaning.h"
#include "iwls.h"
#include "biascorr_beta.h"
#include "offset.h"
#include "stderr.h"


// Non-linear models with fixed effects
// [[Rcpp::export(name = ".bife")]]
Rcpp::List bife(const arma::colvec &y,          // dependant variable
                const arma::mat &X,             // regressor matrix
                const arma::colvec &id,         // column vector with id's
                const arma::colvec &beta_start, // starting values for beta
                const unsigned int model,       // 0 = logit, 1 = probit, 2 = poisson (%MISSING%)
                const unsigned int bias_corr,   // 0 = no, 1 = analytically, 2 = jackknife
                const unsigned int iter_max1,   // iter_max for demeaning 
                const double tolerance1,        // tolerance for demeaning 
                const unsigned int iter_max2,   // iter_max for glm offset
                const double tolerance2) {      // tolerance for glm offset    
 
 
 /* --------------------------------------------------------------------------------------------------------
 *  Prepare data
 * 
 * --------------------------------------------------------------------------------------------------------
 */ 
 
 const arma::uword D = X.n_cols;
 const arma::uword NT = y.n_rows;
 const arma::colvec id_unique = arma::unique(id);
 const arma::uword N = id_unique.n_rows;
 
 arma::colvec T_vector(N);
 arma::colvec y_mean(N);
 arma::mat X_mean(N, D);
 arma::colvec t(NT);
 arma::uvec index_id(NT);
 arma::uword start = 0;
 for (arma::uword i = 0 ; i < N ; ++i) {
  
  arma::uword end = start;
  while (id(end + 1) == id_unique(i)) {
  
   end++;
   if (end == NT - 1) {
    
    break;
   }
  }
  
  const arma::colvec y_i = y.subvec(start, end);
  T_vector(i) = y_i.n_rows;
  y_mean(i) = arma::mean(y_i);
  X_mean.row(i) = arma::mean(X.rows(start, end));
  t.subvec(start, end) = arma::linspace(1, T_vector(i), T_vector(i));
  index_id.subvec(start, end).fill(i);
  
  start = end + 1;
 }
 
 
 /* --------------------------------------------------------------------------------------------------------
 *  Check for linear dependancy
 * 
 * --------------------------------------------------------------------------------------------------------
 */  
 
 if (arma::rank(X) < D) {
  
  Rprintf("\nWarning! The formula contains linear dependant regressor(s).");
  return 0;
 }
 
 
 /* --------------------------------------------------------------------------------------------------------
 *  Compute starting values for alpha
 * 
 * --------------------------------------------------------------------------------------------------------
 */ 
 
 arma::colvec beta = beta_start;
 arma::colvec alpha(N);
 switch (model) {
 
 case 0 : 
  alpha = logistic_inv(y_mean) - X_mean * beta;
  break;
 case 1 :
  alpha = normal_inv(y_mean) - X_mean * beta;
  break;
 }
 
 
 /* --------------------------------------------------------------------------------------------------------
 *  Compute degrees of freedom and standard errors and start demeaning alogrithm
 * 
 * --------------------------------------------------------------------------------------------------------
 */
 
 const unsigned int df = NT - D;
 arma::colvec se_beta(D);
 arma::colvec se_alpha(N);
 arma::mat H_beta(D, D);
 unsigned int iter_demeaning;
 bool conv_demeaning;
 demeaning(beta, alpha, se_beta, se_alpha, H_beta, iter_demeaning, conv_demeaning, y, X, index_id, T_vector, model, iter_max1, tolerance1);
 
 
 /* --------------------------------------------------------------------------------------------------------
 *  Bias Correction 
 *
 *  - No bias correction                                   - bias_corr = 0
 *  - Analytical bias correction - Hahn and Newey (2004)   - bias_corr = 1
 *  - Jackknife bias correction - Hahn and Newey (2004)    - bias_corr = 2 
 * --------------------------------------------------------------------------------------------------------
 */ 
 
 arma::colvec beta_tilde(D);
 switch (bias_corr) {
 
 case 0 : 
  return Rcpp::List::create(Rcpp::Named("par") = Rcpp::List::create(Rcpp::Named("beta") = beta,
                                                                    Rcpp::Named("alpha") = alpha,
                                                                    Rcpp::Named("se.beta") = se_beta,
                                                                    Rcpp::Named("se.alpha") = se_alpha,
                                                                    Rcpp::Named("beta.vcov") = H_beta,
                                                                    Rcpp::Named("avg.alpha") = arma::mean(alpha)),
                            Rcpp::Named("logl.info") = Rcpp::List::create(Rcpp::Named("nobs") = NT,
                                                                          Rcpp::Named("df") = df,
                                                                          Rcpp::Named("loglik") = loglikelihood(y, probability(alpha(index_id) + X * beta, model)),
                                                                          Rcpp::Named("events") = arma::accu(y), 
                                                                          Rcpp::Named("iter.demeaning") = iter_demeaning,
                                                                          Rcpp::Named("conv.demeaning") = conv_demeaning), 
                            Rcpp::Named("model.info") = Rcpp::List::create(Rcpp::Named("used.ids") = id_unique,
                                                                           Rcpp::Named("beta.start") = beta_start,
                                                                           Rcpp::Named("y") = y,
                                                                           Rcpp::Named("X") = X,
                                                                           Rcpp::Named("id") = id,
                                                                           Rcpp::Named("t") = t,
                                                                           Rcpp::Named("index.id") = index_id,
                                                                           Rcpp::Named("T.vector") = T_vector));
 case 1 : 
  beta_tilde = beta - biascorr_beta(beta, alpha, y, X, index_id, T_vector, model);
  break;
 case 2 :
  const arma::uword T = t.max();
  const arma::colvec T_jack = T_vector - 1;
  arma::colvec beta_jack = beta_start;
  arma::colvec alpha_jack(N);
  arma::colvec se_beta_jack(D);
  arma::colvec se_alpha_jack(N);
  arma::mat H_jack(D, D);
  switch (model) {
  
  case 0 : 
   alpha_jack = logistic_inv(y_mean) - X_mean * beta_jack;
   break;
  case 1 :
   alpha_jack = normal_inv(y_mean) - X_mean * beta_jack;
   break;
  }
  unsigned int iter_jack;
  bool conv_jack;
  
  arma::mat beta_mat(D, T);
  for (arma::uword i = 0 ; i < T ; ++i) {
   
   Rcpp::checkUserInterrupt();
   
   const arma::uvec index_jack = arma::find(t != i + 1);
   const arma::colvec y_jack = y(index_jack);
   const arma::mat X_jack = X.rows(index_jack);
   const arma::uvec index_id_jack = index_id(index_jack);
   
   demeaning(beta_jack, alpha_jack, se_beta_jack, se_alpha_jack, H_jack, iter_jack, conv_jack, y_jack, X_jack, index_id_jack, T_jack, model, iter_max1, tolerance1);
   
   beta_mat.col(i) = beta_jack;
  }
  const double mean_T = arma::mean(T_vector);
  beta_tilde = mean_T * beta - (mean_T - 1) * arma::sum(beta_mat, 1) / mean_T;
  break;
 }
 
 /* --------------------------------------------------------------------------------------------------------
 *  Compute alpha_tilde after bias correction
 * 
 * --------------------------------------------------------------------------------------------------------
 */ 
 
 unsigned int iter_offset;
 bool conv_offset;
 const arma::colvec alpha_tilde = glm_offset(iter_offset, conv_offset, beta_tilde, y, X, index_id, T_vector, model, iter_max2, tolerance2);
 
 
 /* --------------------------------------------------------------------------------------------------------
 *  Compute standard errors after bias correction
 * 
 * --------------------------------------------------------------------------------------------------------
 */
 
 arma::colvec se_beta_tilde(D);
 arma::colvec se_alpha_tilde(N);
 arma::mat H_beta_tilde(D, D);
 standard_errors(se_beta_tilde, se_alpha_tilde, H_beta_tilde, beta_tilde, alpha_tilde, y, X, index_id, T_vector, model);
 
 
 return Rcpp::List::create(Rcpp::Named("par") = Rcpp::List::create(Rcpp::Named("beta") = beta,
                                                                   Rcpp::Named("alpha") = alpha,
                                                                   Rcpp::Named("se.beta") = se_beta,
                                                                   Rcpp::Named("se.alpha") = se_alpha,
                                                                   Rcpp::Named("beta.vcov") = H_beta,
                                                                   Rcpp::Named("avg.alpha") = arma::mean(alpha)),
                           Rcpp::Named("par.corr") = Rcpp::List::create(Rcpp::Named("beta") = beta_tilde,
                                                                        Rcpp::Named("alpha") = alpha_tilde,
                                                                        Rcpp::Named("se.beta") = se_beta_tilde,
                                                                        Rcpp::Named("se.alpha") = se_alpha_tilde,
                                                                        Rcpp::Named("beta.vcov") = H_beta_tilde,
                                                                        Rcpp::Named("avg.alpha") = arma::mean(alpha_tilde)),
                           Rcpp::Named("logl.info") = Rcpp::List::create(Rcpp::Named("nobs") = NT,
                                                                         Rcpp::Named("df") = df,
                                                                         Rcpp::Named("loglik") = loglikelihood(y, probability(alpha(index_id) + X * beta, model)),
                                                                         Rcpp::Named("events") = arma::accu(y), 
                                                                         Rcpp::Named("iter.demeaning") = iter_demeaning,
                                                                         Rcpp::Named("conv.demeaning") = conv_demeaning, 
                                                                         Rcpp::Named("loglik.corr") = loglikelihood(y, probability(alpha_tilde(index_id) + X * beta_tilde, model)),
                                                                         Rcpp::Named("iter.offset") = iter_offset,
                                                                         Rcpp::Named("conv.offset") = conv_offset), 
                           Rcpp::Named("model.info") = Rcpp::List::create(Rcpp::Named("used.ids") = id_unique,
                                                                          Rcpp::Named("beta.start") = beta_start,
                                                                          Rcpp::Named("y") = y,
                                                                          Rcpp::Named("X") = X,
                                                                          Rcpp::Named("id") = id,
                                                                          Rcpp::Named("t") = t,
                                                                          Rcpp::Named("index.id") = index_id,
                                                                          Rcpp::Named("T.vector") = T_vector));
}
