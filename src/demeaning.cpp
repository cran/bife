#include "demeaning.h"
#include "iwls.h"


void demeaning(arma::colvec &beta,
               arma::colvec &alpha,
               arma::colvec &se_beta,
               arma::colvec &se_alpha,
               arma::mat &H,
               unsigned int &iter,
               bool &conv,
               const arma::colvec &y,
               const arma::mat &X,
               const arma::uvec &index_id,
               const arma::colvec &T_vector,
               const unsigned int model,
               const unsigned int iter_max,
               const double tolerance) {
 
 
 const arma::uword D = beta.n_rows;
 const arma::uword N = alpha.n_rows;
 const arma::uword NT = y.n_rows;
 
 const arma::colvec step_alpha = {0.0, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0};
 const arma::uword S = step_alpha.n_rows;
 
 arma::colvec swt2(N);
 arma::mat zX(N, D);
 arma::mat X_bar(NT, D);
 
 conv = 0;
 double ll = loglikelihood(y, probability(alpha(index_id) + X * beta, model));
 for (iter = 1 ; iter <= iter_max ; ++iter) {
  
  Rcpp::checkUserInterrupt();
  
  const arma::colvec beta_old = beta;
  const arma::colvec alpha_old = alpha;
  
  const arma::colvec mu = alpha(index_id) + X * beta;
  const arma::colvec w_tilde2 = weight(y, mu, model);
  const arma::colvec w_tilde = arma::sqrt(w_tilde2);
  const arma::mat X_tilde = X.each_col() % w_tilde;
  const arma::colvec y_tilde = gamma(y, mu, model) / w_tilde;
  
  arma::colvec zy(N);
  
  arma::uword start = 0;
  for (arma::uword i = 0 ; i < N ; ++i) {
   
   const arma::uword end = T_vector(i) + start - 1;
   const arma::colvec w_temp = w_tilde.subvec(start, end);
   
   swt2(i) = arma::accu(w_tilde2.subvec(start, end));
   
   zy(i) = arma::accu(w_temp % y_tilde.subvec(start, end)) / swt2(i);
   zX.row(i) = arma::sum(arma::repmat(w_temp, 1, D) % X_tilde.rows(start, end)) / swt2(i);
   
   start = end + 1;
  }
  
  const arma::colvec y_bar = y_tilde - w_tilde % zy(index_id);
  X_bar = X_tilde - arma::repmat(w_tilde, 1, D) % zX.rows(index_id);
  
  H = (X_bar.t() * X_bar).i();
  const arma::colvec beta_diff = H * X_bar.t() * y_bar;
  const arma::colvec alpha_diff = zy - zX * beta_diff;
  
  beta = beta_diff + beta_old;
  alpha = alpha_diff + alpha_old;
  
  const double ll_old = ll;
  ll = loglikelihood(y, probability(alpha(index_id) + X * beta, model));
  if (ll <= ll_old) {
   
   arma::colvec ll_vector(S);
   for (arma::uword i = 0 ; i < S ; ++i) {
    
    const arma::colvec alpha_temp = alpha_old + alpha_diff * step_alpha(i);
    ll_vector(i) = loglikelihood(y, probability(alpha_temp(index_id) + X * beta, model));
   }
   
   const arma::uvec index = arma::find(ll_vector == ll_vector.max());
   alpha = alpha_old + alpha_diff * step_alpha(index(0));
   
   double step_beta = 2.0;
   do{
    
    step_beta *= 0.5;
    const arma::colvec beta_temp = beta_old + beta_diff * step_beta;
    ll = loglikelihood(y, probability(alpha(index_id) + X * beta_temp, model));
   } while (ll < ll_old && step_beta > 1e-6);
   beta = beta_old + beta_diff * step_beta;
  }
  
  const double criterion = arma::norm(beta - beta_old);
  
  if (criterion <= tolerance) {
   
   conv = 1;
   break;
  }
 }
 
 if (conv == 0) {
  
  Rprintf("\nWarning! Maximum number of iterations reached without convergence. (Demeaning)\n");
 }
 
 se_beta = arma::sqrt(H.diag());
 
 for (arma::uword i = 0 ; i < N ; ++i) {
  
  se_alpha(i) = std::sqrt((1.0 / swt2(i)) + arma::as_scalar(zX.row(i) * H * zX.row(i).t()));
 }
}
