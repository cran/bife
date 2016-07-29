#include "link.h"


/* --------------------------------------------------------------------------------------------------------
 *  Logistic Family
 * --------------------------------------------------------------------------------------------------------
 */ 


// Logistic CDF
arma::colvec logistic_cdf(const arma::colvec &x) {
 
 
 return 1.0 / (1.0 + arma::exp(-x));
}


// Logistic PDF
arma::colvec logistic_pdf(const arma::colvec &x) {
 
 
 const arma::colvec temp = arma::exp(-x);
 
 
 return temp / arma::pow(1.0 + temp, 2); 
}


// Logistic invCDF
arma::colvec logistic_inv(const arma::colvec &x) {
 
 
 return arma::log(x / (1.0 - x));
}


// Logistic PDF 1. derivative
arma::colvec deriv1_logistic_pdf(const arma::colvec &x) {
 
 
 const arma::colvec temp = arma::exp(-x);
 
 
 return (temp % (temp - 1.0)) / arma::pow(1.0 + temp, 3);
}


// Logistic PDF 2. derivative
arma::colvec deriv2_logistic_pdf(const arma::colvec &x) {
 
 
 const arma::colvec temp = arma::exp(-x);
 
 
 return (temp % (arma::pow(temp, 2) - 4.0 * temp + 1.0)) / arma::pow(1.0 + temp, 4);
}


/* --------------------------------------------------------------------------------------------------------
*  Normal Family
* --------------------------------------------------------------------------------------------------------
*/ 


// Normal CDF
arma::colvec normal_cdf(const arma::colvec &x) {
 
 
 const arma::uword n = x.n_rows;
 arma::colvec temp(n);
 for (arma::uword i = 0 ; i < n ; ++i) {
  
  temp(i) = R::pnorm(x(i), 0.0, 1.0, 1, 0);
 }
 
 return temp;
}


// Normal PDF
arma::colvec normal_pdf(const arma::colvec &x) {
 
 
 return (1.0 / std::sqrt(2.0 * arma::datum::pi)) * arma::exp(- 0.5 * arma::pow(x, 2));
}


// Normal invCDF
arma::colvec normal_inv(const arma::colvec &x) {
 
 
 const arma::uword n = x.n_rows;
 arma::colvec temp(n);
 for (arma::uword i = 0 ; i < n ; ++i) {
  
  temp(i) = R::qnorm(x(i), 0.0, 1.0, 1, 0);
 }
 
 return temp;
}
