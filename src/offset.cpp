#include "offset.h"
#include "iwls.h"
#include "link.h"


arma::colvec glm_offset(unsigned int &iter,
                        bool &conv,
                        const arma::colvec &beta,
                        const arma::colvec &y,
                        const arma::mat &X,
                        const arma::uvec &index_id,
                        const arma::colvec &T_vector,
                        const unsigned int model,
                        const unsigned int iter_max,
                        const double tolerance) {
 
 
 const arma::uword N = T_vector.n_rows;
 const arma::uword NT = y.n_rows;
 
 const arma::colvec offset = X * beta;
 
 arma::colvec eta(NT);
 switch (model) {
 
 case 0 : 
  eta = logistic_inv((y + 0.5) / 2.0);
  break;
 case 1 :
  eta = normal_inv((y + 0.5) / 2.0);
  break;
 }
 
 conv = 0;
 arma::colvec b = arma::zeros(N);
 for (iter = 1 ; iter <= iter_max ; ++iter) {
  
  Rcpp::checkUserInterrupt();
  
  const arma::colvec b0 = b;
  
  const arma::colvec w = weight(y, eta, model);
  const arma::colvec y_trans = gamma(y, eta, model) / w;
  const arma::colvec z = eta + y_trans - offset;
  const arma::colvec y_adj = z % w;
  
  arma::uword start = 0;
  for (arma::uword i = 0 ; i < N ; ++i) {
   
   const arma::uword end = T_vector(i) + start - 1;
   b(i) = arma::accu(y_adj.subvec(start, end)) / arma::accu(w.subvec(start, end));
   start = end + 1;
  }
  
  if (arma::any(arma::abs(b - b0) < (tolerance * arma::abs(b0)))) {
   
   conv = 1;
   break;
  }
  
  eta = b(index_id) + offset;
 }
 
 if (conv == 0) {
  
  Rprintf("\nWarning! Maximum number of iterations reached without convergence. (Offset)\n");
 }
 
 
 return b;
}
