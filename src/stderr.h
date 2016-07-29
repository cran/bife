#define ARMA_NO_DEBUG


#include <RcppArmadillo.h>


#ifndef __STDERR__
#define __STDERR__

void standard_errors(arma::colvec &se_beta, arma::colvec &se_alpha, arma::mat &H,
                     const arma::colvec &beta, const arma::colvec &alpha,
                     const arma::colvec &y, const arma::mat &X,
                     const arma::uvec &index_id, const arma::colvec &T_vector,
                     const unsigned int model);

#endif // __STDERR__
