#include "biascorr_apeff.h"
#include "iwls.h"
#include "link.h"


// Analytical bias-correction for average partial effects
arma::colvec biascorr_apeff(const arma::colvec &discrete,
                            const arma::colvec &beta_tilde,
                            const arma::colvec &alpha_tilde,
                            const arma::colvec &y,
                            const arma::mat &X,
                            const arma::uvec &index_id,
                            const arma::colvec &T_vector,
                            const unsigned int model) {
 
 
 const arma::uword D = beta_tilde.n_rows;
 const arma::uword N = alpha_tilde.n_rows;
 const arma::uword NT = y.n_rows;
 
 const arma::colvec mu = X * beta_tilde + alpha_tilde(index_id);
 const arma::colvec v_it = gamma(y, mu, model);
 const arma::colvec v_ita = - weight(y, mu, model);
 const arma::colvec v_it2 = arma::pow(v_it, 2);
 const arma::colvec V2_it = v_it2 + v_ita;
 
 arma::colvec sigma_i2(NT);
 arma::colvec beta_i(NT);
 arma::uword start = 0;
 for (arma::uword i = 0 ; i < N ; ++i) {
  
  const arma::uword end = T_vector(i) + start - 1;
  sigma_i2.subvec(start, end).fill(T_vector(i) / arma::accu(v_it2.subvec(start, end)));
  beta_i.subvec(start, end) = - arma::pow(sigma_i2.subvec(start, end), 2) * arma::accu(v_it.subvec(start, end) % V2_it.subvec(start, end)) / (2.0 * T_vector(i));
  start = end + 1;
 }
 
 const arma::colvec phi_it = sigma_i2 % v_it;
 
 arma::mat m_a(NT, D);
 arma::mat m_aa(NT, D);
 for (arma::uword i = 0 ; i < D ; ++i) {
  
  if (discrete(i) == 1) {
   
   arma::mat X_one = X;
   arma::mat X_zero = X;
   X_one.col(i).ones();
   X_zero.col(i).zeros();
   const arma::colvec mu_one = X_one * beta_tilde + alpha_tilde(index_id);
   const arma::colvec mu_zero = X_zero * beta_tilde + alpha_tilde(index_id);
   
   switch (model) {
   
   case 0 :
    m_a.col(i) = logistic_pdf(mu_one) - logistic_pdf(mu_zero);
    m_aa.col(i) = deriv1_logistic_pdf(mu_one) - deriv1_logistic_pdf(mu_zero);
    break;
   case 1 :
    m_a.col(i) = normal_pdf(mu_one) - normal_pdf(mu_zero);
    m_aa.col(i) = - mu_one % normal_pdf(mu_one) + mu_zero % normal_pdf(mu_zero);
    break;
   }
  } else {
   
   switch (model) {
   
   case 0 :
    m_a.col(i) = beta_tilde(i) * deriv1_logistic_pdf(mu);
    m_aa.col(i) = beta_tilde(i) * deriv2_logistic_pdf(mu);
    break;
   case 1 :
    m_a.col(i) = - beta_tilde(i) * mu % normal_pdf(mu);
    m_aa.col(i) = beta_tilde(i) * normal_pdf(mu) % (arma::pow(mu, 2) - 1.0) ;
    break;
   }
  }
 }
 
 const arma::rowvec delta_hat = arma::mean(m_a.each_col() % (beta_i + phi_it) + m_aa.each_col() % (sigma_i2 / 2.0));
 
 
 return delta_hat.t() / arma::mean(T_vector);
}
