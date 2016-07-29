#define ARMA_NO_DEBUG


#include <RcppArmadillo.h>


#ifndef __BIASCORR_APEFF__
#define __BIASCORR_APEFF__

arma::colvec biascorr_apeff(const arma::colvec &discrete, const arma::colvec &beta_tilde,
                            const arma::colvec &alpha_tilde, const arma::colvec &y,
                            const arma::mat &X, const arma::uvec &index_id,
                            const arma::colvec &T_vector, const unsigned int model);

#endif // __BIASCORR_APEFF__
