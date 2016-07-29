#define ARMA_NO_DEBUG


#include <RcppArmadillo.h>


#ifndef __IWLS__
#define __IWLS__

arma::colvec probability(const arma::colvec &mu, const unsigned int model);
double loglikelihood(const arma::colvec &y, const arma::colvec &p);
arma::colvec gamma(const arma::colvec &y, const arma::colvec &mu,
                   const unsigned int model);
arma::colvec weight(const arma::colvec &y, const arma::colvec &mu,
                    const unsigned int model);

#endif // __IWLS__
