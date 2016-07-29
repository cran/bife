#define ARMA_NO_DEBUG


#include <RcppArmadillo.h>


#ifndef __DEMEANING__
#define __DEMEANING__

void demeaning(arma::colvec &beta, arma::colvec &alpha,
               arma::colvec &se_beta, arma::colvec &se_alpha,
               arma::mat &H, unsigned int &iter, bool &conv,
               const arma::colvec &y, const arma::mat &X,
               const arma::uvec &index_id, const arma::colvec &T_vector,
               const unsigned int model, const unsigned int iter_max,
               const double tolerance);

#endif // __DEMEANING__
