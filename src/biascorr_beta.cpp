#include "biascorr_beta.h"
#include "iwls.h"


// Analytical bias-correction for betas
arma::colvec biascorr_beta(const arma::colvec &beta,
                           const arma::colvec &alpha,
                           const arma::colvec &y,
                           const arma::mat &X,
                           const arma::uvec &index_id,
                           const arma::colvec &T_vector,
                           const unsigned int model) {
 
 
 // Dimensions of arrays
 const arma::uword D = X.n_cols;
 const arma::uword N = alpha.n_rows;
 const arma::uword NT = y.n_rows;
 
 // Compute u_it, v_it, and v_ita
 const arma::colvec mu = X * beta + alpha(index_id);
 const arma::colvec v_it = gamma(y, mu, model);
 const arma::mat u_it = X.each_col() % v_it;
 const arma::colvec v_ita = - weight(y, mu, model);
 
 // Compute v_it2, uv, and V2_it
 const arma::colvec v_it2 = arma::pow(v_it, 2);
 const arma::colvec V2_it = v_it2 + v_ita;
 const arma::mat uv = u_it.each_col() % v_it;
 
 // Compute sum of u_it.each_col() % v_it
 arma::uword start;
 arma::mat sum_uv(NT, D);
 for (arma::uword j = 0 ; j < D ; ++j) {
  
  start = 0;
  for (arma::uword i = 0 ; i < N ; ++i) {
   
   const arma::uword end = T_vector(i) + start - 1;
   sum_uv(arma::span(start, end), j).fill(arma::accu(uv(arma::span(start, end), j)));
   start = end + 1;
  }
 }
 
 // Compute sum of v_it2
 arma::colvec sum_v2(NT);
 start = 0;
 for (arma::uword i = 0 ; i < N ; ++i) {
  
  const arma::uword end = T_vector(i) + start - 1;
  sum_v2.subvec(start, end).fill(arma::accu(v_it2.subvec(start, end)));
  start = end + 1;
 }
 
 // Compute U_it
 const arma::mat temp = sum_uv.each_col() / sum_v2;
 const arma::mat U_it = u_it - v_it % temp.each_col();
 
 // Compute B_bar, H_bar, and b_bar
 const arma::mat H_bar = (U_it.each_col() / T_vector(index_id)).t() * U_it;
 const arma::colvec b_bar = 0.5 * U_it.t() * (V2_it / sum_v2);
 const arma::colvec B_bar = - H_bar.i() * b_bar;
 
 
 return B_bar / arma::mean(T_vector);
}
