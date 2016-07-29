#include "stderr.h"
#include "iwls.h"


void standard_errors(arma::colvec &se_beta,
                     arma::colvec &se_alpha,
                     arma::mat &H,
                     const arma::colvec &beta,
                     const arma::colvec &alpha,
                     const arma::colvec &y,
                     const arma::mat &X,
                     const arma::uvec &index_id,
                     const arma::colvec &T_vector,
                     const unsigned int model) {
 
 
 const arma::uword N = alpha.n_rows;
 const arma::uword D = beta.n_rows;
 const arma::colvec w_tilde2 = weight(y, alpha(index_id) + X * beta, model);
 const arma::colvec w_tilde = arma::sqrt(w_tilde2);
 const arma::mat X_tilde = X.each_col() % w_tilde;
 
 arma::mat zX(N, D);
 arma::colvec swt2(N);
 arma::uword start = 0;
 for (arma::uword i = 0 ; i < N ; ++i) {
  
  const arma::uword end = T_vector(i) + start - 1;
  
  swt2(i) = arma::accu(w_tilde2.subvec(start, end));
  
  zX.row(i) = arma::sum(arma::repmat(w_tilde.subvec(start, end), 1, D) % X_tilde.rows(start, end)) / swt2(i);
  
  start = end + 1;
 }
 
 const arma::mat X_bar = X_tilde - arma::repmat(w_tilde, 1, D) % zX.rows(index_id);
 
 H = (X_bar.t() * X_bar).i();
 se_beta = arma::sqrt(H.diag());
 
 for (arma::uword i = 0 ; i < N ; ++i) {
  
  se_alpha(i) = std::sqrt((1.0 / swt2(i)) + arma::as_scalar(zX.row(i) * H * zX.row(i).t()));
 }
}
