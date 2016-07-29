#include "iwls.h"
#include "link.h"


// Probability for likelihood function
arma::colvec probability(const arma::colvec &mu,
                         const unsigned int model) {
 
 
 const arma::uword n = mu.n_rows;
 
 arma::colvec p(n);
 switch (model) {
 
 case 0 :
  p = logistic_cdf(mu);
  break;
 case 1 :
  p = normal_cdf(mu);
  break;
 }
 
 p(arma::find(p == 0)).fill(1e-15);
 p(arma::find(p == 1)).fill(1.0 - 1e-15);
 
 
 return p;
}


// Loglikelihood function
double loglikelihood(const arma::colvec &y,
                     const arma::colvec &p) {
 
 
 const arma::uvec index = arma::find(y == 0);
 arma::colvec P = p;
 P(index) = 1.0 - p(index);
 
 
 return arma::accu(arma::log(P));
}


// "gamma" for IWLS / Gradient
arma::colvec gamma(const arma::colvec &y,
                   const arma::colvec &mu,
                   const unsigned int model) {
 
 
 const arma::uword n = y.n_rows;
 
 arma::colvec g(n);
 switch (model) {
 
 case 0 :
  g = y - logistic_cdf(mu);
  break;
 case 1 :
  const arma::colvec p = probability(mu, model);
  g = (normal_pdf(mu) % (y - normal_cdf(mu))) / (p % (1.0 - p));
  break;
 }
 
 
 return g;
}


// W of IWLS / Hessian
arma::colvec weight(const arma::colvec &y,
                    const arma::colvec &mu,
                    const unsigned int model) {
 
 
 const arma::uword n = y.n_rows;
 
 arma::colvec w(n);
 switch (model) {
 
 case 0 :
  w = logistic_pdf(mu);
  w(arma::find(w == 0)).fill(1e-15);
  break;
 case 1 :
  const arma::colvec p = probability(mu, model);
  const arma::colvec h_it = normal_pdf(mu) / (p % (1.0 - p));
  const arma::colvec dh_it = mu % h_it + (arma::pow(normal_pdf(mu), 2) % (1.0 - 2.0 * normal_cdf(mu))) / arma::pow(p % (1.0 - p), 2);
  w = dh_it % (y - normal_cdf(mu)) + normal_pdf(mu) % h_it;
  break;
 }
 
 
 return w;
}
