#define ARMA_NO_DEBUG


#include <RcppArmadillo.h>


#ifndef __BIASCORR_BETA__
#define __BIASCORR_BETA__

arma::colvec biascorr_beta(const arma::colvec &beta, const arma::colvec &alpha,
                           const arma::colvec &y, const arma::mat &X,
                           const arma::uvec &index_id, const arma::colvec &T_vector,
                           const unsigned int model);

#endif // __BIASCORR_BETA__
