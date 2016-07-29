#define ARMA_NO_DEBUG


#include <RcppArmadillo.h>


#ifndef __OFFSET__
#define __OFFSET__

arma::colvec glm_offset(unsigned int &iter, bool &conv, const arma::colvec &beta,
                        const arma::colvec &y, const arma::mat &X,
                        const arma::uvec &index_id, const arma::colvec &T_vector,
                        const unsigned int model, const unsigned int iter_max,
                        const double tolerance);

#endif // __OFFSET__
