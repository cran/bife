#define ARMA_NO_DEBUG


#include <RcppArmadillo.h>


#ifndef __LINK__
#define __LINK__

arma::colvec logistic_cdf(const arma::colvec &x);
arma::colvec logistic_pdf(const arma::colvec &x);
arma::colvec logistic_inv(const arma::colvec &x);
arma::colvec deriv1_logistic_pdf(const arma::colvec &x);
arma::colvec deriv2_logistic_pdf(const arma::colvec &x);
arma::colvec normal_cdf(const arma::colvec &x);
arma::colvec normal_pdf(const arma::colvec &x);
arma::colvec normal_inv(const arma::colvec &x);

#endif // __LINK__
