#define ARMA_NO_DEBUG


#include <RcppArmadillo.h>



/* --- CDF, PDF, and ICDF of logistic and normal distribution
 * 
 * List of functions:
 * 
 *    name            description
 *    -----------------------------------------------------------------------------
 *    PDF           - Corresponding probability density function
 *    CDF           - Corresponding cumulative distribution function
 *    ICDF          - Inverse of the corresponding cumulative distribution function
 *    DPDF_alpha    - Corresponding first derivative of PDF() wrt to alpha
 *    DDPDF_alpha   - Corresponding second derivative of PDF() wrt to alpha
 *    
 */


arma::colvec PDF(const arma::colvec &x,
                 const unsigned int model) {
  
  
  const arma::uword n = x.n_rows;
  
  arma::colvec pdf(n);
  switch (model) {
  
    case 0 :
      {
        for (arma::uword i = 0 ; i < n ; ++i) pdf(i) = ::Rf_dlogis(x(i), 0.0, 1.0, 0);
        break;
      }
    case 1 :
      {
        for (arma::uword i = 0 ; i < n ; ++i) pdf(i) = ::Rf_dnorm4(x(i), 0.0, 1.0, 0);
        break;
      }
  }
  

  return pdf;
}


arma::colvec CDF(const arma::colvec &x,
                 const unsigned int model) {
  
  
  const arma::uword n = x.n_rows;
  
  arma::colvec cdf(n);
  switch (model) {
  
    case 0 :
      {
        for (arma::uword i = 0 ; i < n ; ++i) cdf(i) = ::Rf_plogis(x(i), 0.0, 1.0, 1, 0);
        break;
      }
    case 1 :
      {
        for (arma::uword i = 0 ; i < n ; ++i) cdf(i) = ::Rf_pnorm5(x(i), 0.0, 1.0, 1, 0);
        break;  
      }
  }
  
  
  return cdf;
}


arma::colvec ICDF(const arma::colvec &x,
                  const unsigned int model) {
  
  
  const arma::uword n = x.n_rows;
  
  arma::colvec icdf(n);
  switch (model) {
  
    case 0 :
      {
        for (arma::uword i = 0 ; i < n ; ++i) icdf(i) = ::Rf_qlogis(x(i), 0.0, 1.0, 1, 0);
        break;
      }
    case 1 :
      {
        for (arma::uword i = 0 ; i < n ; ++i) icdf(i) = ::Rf_qnorm5(x(i), 0.0, 1.0, 1, 0);
        break;
      }
  }
  
  
  return icdf;
}


arma::colvec DPDF_alpha(const arma::colvec &x,
                        const unsigned int model) {
  
  
  const arma::uword n = x.n_rows;
  
  arma::colvec dpdf(n);
  switch (model) {
  
    case 0 :
      {
        dpdf = PDF(x, model) % (1.0 - 2.0 * CDF(x, model));
        break;
      }
    case 1 :
      {
        dpdf = - x % PDF(x, model);
        break;
      }
  }
  
  
  return dpdf;
}


arma::colvec DDPDF_alpha(const arma::colvec &x,
                         const unsigned int model) {
  
  
  const arma::uword n = x.n_rows;
  
  arma::colvec ddpdf(n);
  switch (model) {
  
    case 0 :
      {
        const arma::colvec p = PDF(x, model);
        ddpdf = p % (arma::pow(1.0 - 2.0 * CDF(x, model), 2) - 2.0 * p);
        break;
      }
    case 1 :
      {
        ddpdf = PDF(x, model) % (arma::pow(x, 2) - 1.0);
        break;
      }
  }
  
  
  return ddpdf;
}


/* --- Functions for IWLS algorithm
 * 
 * List of functions:
 * 
 *    name            description
 *    -----------------------------------------------------------------------------
 *    log_likelihood  - Computes the value of the log-likelihood 
 *    gamma           - Computes 'gamma' of the IWLS algorithm (gradient = X' * gamma)
 *    weights         - Computes the 'weights' of the IWLS algorithm (Hessian = - X' * W * X)
 *    
 */


double log_likelihood(const arma::colvec &z,
                      const unsigned int model) {
  
  
  arma::colvec cdf = CDF(z, model);
  cdf(arma::find(cdf == 0.0)).fill(1.0e-16);
  
  
  return arma::accu(arma::log(cdf));
}


// Computes 'gamma' such that the gradient equals X' * gamma
arma::colvec gamma(const arma::colvec &y,
                   const arma::colvec &mu,
                   const unsigned int model) {
  
  
  const arma::uword n = y.n_rows;
  
  arma::colvec g(n);
  switch (model) {
  
    case 0 :
      {
        g = y - CDF(mu, model);
        break;
      }
    case 1 :
      {
        const arma::colvec q = 2.0 * y - 1.0;
        const arma::colvec z = q % mu;
        arma::colvec Phi = CDF(z, model);
        Phi(arma::find(Phi == 0.0)).fill(1.0e-16);
        
        g = (PDF(z, model) % q) / Phi;
        break;
      }
  }
  
  
  return g;
}


// Computes the 'weights' of the diagonal matrix that is needed to compute the Hessian (- X' * W * X)
arma::colvec weights(const arma::colvec &y,
                     const arma::colvec &mu,
                     const unsigned int model) {
  
  
  const arma::uword n = y.n_rows;
  
  arma::colvec w(n);
  switch (model) {
  
    case 0 :
      {
        arma::colvec pdf = PDF(mu, model);
        pdf(arma::find(pdf == 0.0)).fill(1.0e-16);
        
        w = pdf;
        break;
      }
    case 1 :
      {
        const arma::colvec z = (2.0 * y - 1.0) % mu;
        const arma::colvec phi = PDF(z, model);
        const arma::colvec Phi = CDF(z, model);
        w = phi % (phi + z % Phi) / arma::pow(Phi, 2);
        
        w(arma::find(w == 0.0)).fill(1.0e-16);
        break;
      }
    }
  
  
  return w;
}


/* --- Function demeaning()
 * 
 * Input arguments:
 * 
 *    name          dim       description
 *    -----------------------------------------------------------------------------
 *    %Missing%
 * 
 */


void demeaning(arma::colvec &beta,
               arma::colvec &alpha,
               arma::colvec &se_beta,
               arma::colvec &se_alpha,
               arma::mat &H,
               unsigned int &iter,
               bool &conv,
               const arma::colvec &y,
               const arma::mat &X,
               const arma::uvec &index_id,
               const arma::colvec &Ti_vector,
               const unsigned int model,
               const unsigned int iter_max,
               const double tolerance) {
  
  
  // Preliminary
  const arma::uword d = beta.n_rows;
  const arma::uword n = alpha.n_rows;
  const arma::uword s = 8;
  
  // Step size adjustment for 'alpha'
  arma::rowvec step_alpha(s);
  step_alpha << 0.0 << 0.0625 << 0.125 << 0.25 << 0.5 << 1.0 << 2.0 << 4.0 << arma::endr;
  
  // Compute 'q' for the loglikelihood
  const arma::colvec q = 2.0 * y - 1.0;
  
  // 'Demeaning' algorithm
  arma::colvec swt_sq(n);
  arma::mat XX(n, d);
  conv = 0;
  double ll = log_likelihood(q % (alpha(index_id) + X * beta), model);
  for (iter = 1 ; iter <= iter_max ; ++iter) {
    
    Rcpp::checkUserInterrupt();
    
    const arma::colvec beta_old = beta;
    const arma::colvec alpha_old = alpha;
    const arma::colvec mu = alpha(index_id) + X * beta;
    const arma::colvec w_tilde_sq = weights(y, mu, model);
    const arma::colvec w_tilde = arma::sqrt(w_tilde_sq);
    const arma::mat X_tilde = X.each_col() % w_tilde;
    const arma::colvec y_tilde = gamma(y, mu, model) / w_tilde;
    
    arma::colvec yy(n);
    arma::uword start = 0;
    for (arma::uword i = 0 ; i < n ; ++i) {
      
      const arma::uword end = Ti_vector(i) + start - 1;
      const arma::colvec wi_tilde = w_tilde.subvec(start, end);
      
      swt_sq(i) = arma::accu(w_tilde_sq.subvec(start, end));
      yy(i) = arma::accu(wi_tilde % y_tilde.subvec(start, end)) / swt_sq(i);
      XX.row(i) = arma::sum(arma::repmat(wi_tilde, 1, d) % X_tilde.rows(start, end)) / swt_sq(i);
      
      start = end + 1;
    }
    
    const arma::colvec y_bar = y_tilde - w_tilde % yy(index_id);
    const arma::mat X_bar = X_tilde - arma::repmat(w_tilde, 1, d) % XX.rows(index_id);
    
    H = (X_bar.t() * X_bar).i();
    const arma::colvec beta_diff = H * X_bar.t() * y_bar;
    const arma::colvec alpha_diff = yy - XX * beta_diff;
    beta = beta_diff + beta_old;
    alpha = alpha_diff + alpha_old;
    
    const double ll_old = ll;
    ll = log_likelihood(q % (alpha(index_id) + X * beta), model);
    if (ll <= ll_old) {
      
      double step_alphai;
      for (arma::uword i = 0 ; i < s ; ++i) {
        
        const arma::colvec alpha_temp = alpha_old + alpha_diff * step_alpha(i);
        const double ll_temp = log_likelihood(q % (alpha_temp(index_id) + X * beta), model);
        if (ll_temp >= ll) {
          
          ll = ll_temp;
          step_alphai = step_alpha(i);
        }
      }
      
      alpha = alpha_old + alpha_diff * step_alphai;
   
      double step_beta = 2.0;
      do{
        
        step_beta *= 0.5;
        const arma::colvec beta_temp = beta_old + beta_diff * step_beta;
        ll = log_likelihood(q % (alpha(index_id) + X * beta_temp), model);
      } while (ll < ll_old && step_beta > 1.0e-06);
      
      beta = beta_old + beta_diff * step_beta;
    }
    
    if (arma::norm(beta - beta_old) <= tolerance) {
      
      conv = 1;
      break;
    }
  }
  
  if (conv == 0) {
    
    Rcpp::Rcout << "\nWarning! Maximum number of iterations reached without convergence. (Demeaning)\n" << std::endl;
  }
  
  se_beta = arma::sqrt(H.diag());
  for (arma::uword i = 0 ; i < n ; ++i) se_alpha(i) = std::sqrt((1.0 / swt_sq(i)) + arma::as_scalar(XX.row(i) * H * XX.row(i).t()));
}


/* --- Function offset()
 * 
 * Input arguments:
 * 
 *    name          dim       description
 *    -----------------------------------------------------------------------------
 *    %Missing%
 * 
 */


arma::colvec offset(unsigned int &iter,
                    bool &conv,
                    const arma::colvec &beta,
                    const arma::colvec &y,
                    const arma::mat &X,
                    const arma::uvec &index_id,
                    const arma::colvec &Ti_vector,
                    const unsigned int model,
                    const unsigned int iter_max,
                    const double tolerance) {
  
  
  const arma::uword n = Ti_vector.n_rows;
  const arma::colvec mu = X * beta;

  conv = 0;
  arma::colvec eta = ICDF((y + 0.5) / 2.0, model);
  arma::colvec alpha_adj = arma::zeros(n);
  for (iter = 1 ; iter <= iter_max ; ++iter) {
    
    Rcpp::checkUserInterrupt();
    
    const arma::colvec alpha_adj_old = alpha_adj;
    const arma::colvec w = weights(y, eta, model);
    const arma::colvec y_adj = (eta + gamma(y, eta, model) / w - mu) % w;
    
    arma::uword start = 0;
    for (arma::uword i = 0 ; i < n ; ++i) {
    
      const arma::uword end = Ti_vector(i) + start - 1;
      
      alpha_adj(i) = arma::accu(y_adj.subvec(start, end)) / arma::accu(w.subvec(start, end));
      
      start = end + 1;
    }
    
    if (arma::any(arma::abs(alpha_adj - alpha_adj_old) < (tolerance * arma::abs(alpha_adj_old)))) {
      
      conv = 1;
      break;
    }
    
    eta = alpha_adj(index_id) + mu;
  }
  
  if (conv == 0) {
    
    Rprintf("\nWarning! Maximum number of iterations reached without convergence. (Offset)\n");
  }
  
  
  return alpha_adj;
}


/* --- Function standard_errors()
 * 
 * Input arguments:
 * 
 *    name          dim       description
 *    -----------------------------------------------------------------------------
 *    %Missing%
 * 
 */


void standard_errors(arma::colvec &se_beta,
                     arma::colvec &se_alpha,
                     arma::mat &H,
                     const arma::colvec &beta,
                     const arma::colvec &alpha,
                     const arma::colvec &y,
                     const arma::mat &X,
                     const arma::uvec &index_id,
                     const arma::colvec &Ti_vector,
                     const unsigned int model) {
  
  
  const arma::uword n = alpha.n_rows;
  const arma::uword d = beta.n_rows;
  const arma::colvec w_tilde_sq = weights(y, alpha(index_id) + X * beta, model);
  const arma::colvec w_tilde = arma::sqrt(w_tilde_sq);
  const arma::mat X_tilde = X.each_col() % w_tilde;
  
  arma::mat XX(n, d);
  arma::colvec swt_sq(n);
  arma::uword start = 0;
  for (arma::uword i = 0 ; i < n ; ++i) {
    
    const arma::uword end = Ti_vector(i) + start - 1;
    
    swt_sq(i) = arma::accu(w_tilde_sq.subvec(start, end));
    XX.row(i) = arma::sum(arma::repmat(w_tilde.subvec(start, end), 1, d) % X_tilde.rows(start, end)) / swt_sq(i);
    
    start = end + 1;
  }
  
  const arma::mat X_bar = X_tilde - arma::repmat(w_tilde, 1, d) % XX.rows(index_id);
  
  H = (X_bar.t() * X_bar).i();
  se_beta = arma::sqrt(H.diag());
  
  for (arma::uword i = 0 ; i < n ; ++i) se_alpha(i) = std::sqrt((1.0 / swt_sq(i)) + arma::as_scalar(XX.row(i) * H * XX.row(i).t()));
}


/* --- Functions to compute average partial effects
 * 
 * List of functions:
 * 
 *    name            description
 *    -----------------------------------------------------------------------------
 *    
 *    avg_peff        - Computes average partial effects 
 *    avg_peff_corr   - Computes corrected average partial effects
 * 
 * Input arguments:
 * 
 *    name          dim       description
 *    -----------------------------------------------------------------------------
 *    %Missing%
 * 
 */


// Computes average partial effects
// [[Rcpp::export(name = ".avg_peff")]]
arma::colvec avg_peff(const arma::colvec &discrete,
                      const arma::colvec &beta,
                      const arma::colvec &alpha,
                      const arma::mat &X,
                      const arma::uvec &index_id,
                      const unsigned int model) {
  
  
  const arma::uword d = X.n_cols;
  const arma::uword nt = X.n_rows;
  
  const arma::colvec mu = X * beta + alpha(index_id);
  arma::mat ind_peff(nt, d);
  for (arma::uword i = 0 ; i < d ; ++i) {
    
    if (discrete(i) == 1) {
      
      arma::mat X_inc = X;
      X_inc.col(i) += 1.0;
      ind_peff.col(i) = CDF(X_inc * beta + alpha(index_id), model) - CDF(mu, model);
    } else {
      
      ind_peff.col(i) = beta(i) * PDF(mu, model);
    }
  }
  
  
  return arma::mean(ind_peff).t();
}


// Corrected average partial effects
// [[Rcpp::export(name = ".avg_peff_corr")]]
arma::mat avg_peff_corr(const arma::colvec &discrete,
                        const arma::colvec &beta,
                        const arma::colvec &alpha,
                        const arma::colvec &beta_start,
                        const arma::colvec &alpha_start,
                        const arma::colvec &beta_tilde,
                        const arma::colvec &alpha_tilde,
                        const arma::colvec &y,
                        const arma::mat &X,
                        const arma::colvec &time,
                        const arma::uvec &index_id,
                        const arma::colvec &Ti_vector,
                        const unsigned int model,
                        const unsigned int bias_corr,
                        const unsigned int iter_max1,  
                        const double tolerance1,
                        const unsigned int iter_max2,  
                        const double tolerance2) {    
  
  
  const arma::uword n = alpha.n_rows;
  const arma::uword d = beta.n_rows;
  const arma::uword nt = y.n_rows;
  
  arma::mat Avg_peff(d, 2);
  Avg_peff.col(0) = avg_peff(discrete, beta, alpha, X, index_id, model);
  switch (bias_corr) {
  
    case 0 :
      {
        // 0 - Semi bias correction
        Avg_peff.col(1) = avg_peff(discrete, beta_tilde, alpha_tilde, X, index_id, model);
        break;
      }
    case 1 :
      {
        // 1 - Analytical bias correction
        const arma::colvec mu = X * beta_tilde + alpha_tilde(index_id);
        const arma::colvec vit = gamma(y, mu, model);
        const arma::colvec vit_alpha = - weights(y, mu, model);
        const arma::colvec vit_sq = arma::pow(vit, 2);
        const arma::colvec V2it = vit_sq + vit_alpha;
        
        arma::colvec sigmai_sq(nt);
        arma::colvec betai(nt);
        arma::uword start = 0;
        for (arma::uword i = 0 ; i < n ; ++i) {
          
          const arma::uword end = Ti_vector(i) + start - 1;
          
          sigmai_sq.subvec(start, end).fill(Ti_vector(i) / arma::accu(vit_sq.subvec(start, end)));
          betai.subvec(start, end) = - arma::pow(sigmai_sq.subvec(start, end), 2) * arma::accu(vit.subvec(start, end) % V2it.subvec(start, end)) / (2.0 * Ti_vector(i));
          
          start = end + 1;
        }
        
        const arma::colvec phiit = sigmai_sq % vit;
        
        arma::mat m_a(nt, d);
        arma::mat m_aa(nt, d);
        for (arma::uword i = 0 ; i < d ; ++i) {
          
          if (discrete(i) == 1) {
            
            arma::mat X_inc = X;
            X_inc.col(i) += 1.0;
            const arma::colvec mu_inc = X_inc * beta_tilde + alpha_tilde(index_id);
            
            m_a.col(i) = PDF(mu_inc, model) - PDF(mu, model);
            m_aa.col(i) = DPDF_alpha(mu_inc, model) - DPDF_alpha(mu, model);
          } else {
            
            m_a.col(i) = beta_tilde(i) * DPDF_alpha(mu, model);
            m_aa.col(i) = beta_tilde(i) * DDPDF_alpha(mu, model);
          }
        }
        
        const arma::colvec delta1 = m_a.t() * (betai + phiit);
        const arma::colvec delta2 = m_aa.t() * (sigmai_sq / 2.0);
        const arma::colvec delta_hat = (delta1 + delta2) / nt;
        
        Avg_peff.col(1) = avg_peff(discrete, beta_tilde, alpha_tilde, X, index_id, model) - delta_hat / arma::mean(Ti_vector);
        break;
      }
  }
  
  
  return Avg_peff;
}


/* --- Function bife()
 * 
 * Input arguments:
 * 
 *    name          dim       description
 *    -----------------------------------------------------------------------------
 *    y           - nt x 1  - dependant variable (binary)
 *    X           - nt x d  - regressor matrix
 *    id          - nt x 1  - individual identifier 
 *    beta_start  - d x 1   - starting values of 'beta' (structural parameters)
 *    model       - scalar  - 0 = logit, 1 = probit
 *    bias_corr   - scalar  - 0 = no, 1 = analytically, 2 = jackknife
 *    iter_max1   - scalar  - max # iterations (Demeaning)
 *    tolerance1  - scalar  - tolerance level (Demeaning)
 *    iter_max2   - scalar  - max # iterations (Offset)
 *    tolerance2  - scalar  - tolerance level (Offset)
 * 
 */


// [[Rcpp::export(name = ".bife")]]
Rcpp::List bife(const arma::colvec &y,
                const arma::mat &X,
                const arma::colvec &id,
                const arma::colvec &beta_start,
                const unsigned int model,
                const unsigned int bias_corr,
                const unsigned int iter_max1,
                const double tolerance1,
                const unsigned int iter_max2,
                const double tolerance2) {  
  
 
  // Prepare data
  const arma::uword d = X.n_cols;
  const arma::uword nt = y.n_rows;
  const arma::colvec id_unique = arma::unique(id);
  const arma::uword n = id_unique.n_rows;
  
  arma::colvec Ti_vector(n);
  arma::colvec y_bar(n);
  arma::mat X_bar(n, d);
  arma::colvec time(nt);
  arma::uvec index_id(nt);
  arma::uword start = 0;
  for (arma::uword i = 0 ; i < n ; ++i) {
    
    arma::uword end = start;
    while (end < nt && id(++end) == id_unique(i));
    const arma::colvec yi = y.subvec(start, --end);
    const arma::uword Ti = yi.n_rows;
    
    Ti_vector(i) = Ti;
    y_bar(i) = arma::mean(yi);
    X_bar.row(i) = arma::mean(X.rows(start, end));
    time.subvec(start, end) = arma::linspace(1, Ti, Ti);
    index_id.subvec(start, end).fill(i);
    
    start = end + 1;
  }
  
  // Check for full column rank
  if (arma::rank(X) < d) {
    
    Rcpp::Rcout << "\nWarning! The formula contains linear dependant regressor(s)." << std::endl;
    return 0;
  }
  
  // Compute starting values for 'alpha'
  arma::colvec beta = beta_start;
  const arma::colvec alpha_start = ICDF(y_bar, model) - X_bar * beta;
  arma::colvec alpha = alpha_start;
  
  // Start 'demeaning' to obtain coefficients and standard errors
  arma::colvec se_beta(d);
  arma::colvec se_alpha(n);
  arma::mat H_beta(d, d);
  unsigned int iter_demeaning;
  bool conv_demeaning;
  demeaning(beta, alpha, se_beta, se_alpha, H_beta, iter_demeaning, conv_demeaning, y, X, index_id, Ti_vector, model, iter_max1, tolerance1);
  
  // Apply bias correction if bias_corr > 0
  arma::colvec beta_tilde(d);
  switch (bias_corr) {
  
    case 0 :
      {
        // 0 - No bias correction
        
        // Return "0 - No bias correction"
        return Rcpp::List::create(Rcpp::Named("par") = Rcpp::List::create(Rcpp::Named("beta") = beta,
                                                                          Rcpp::Named("alpha") = alpha,
                                                                          Rcpp::Named("se_beta") = se_beta,
                                                                          Rcpp::Named("se_alpha") = se_alpha,
                                                                          Rcpp::Named("beta_vcov") = H_beta,
                                                                          Rcpp::Named("avg_alpha") = arma::mean(alpha)),
                                  Rcpp::Named("logl_info") = Rcpp::List::create(Rcpp::Named("loglik") = log_likelihood((2.0 * y - 1.0) % (X * beta + alpha(index_id)), model),
                                                                                Rcpp::Named("iter_demeaning") = iter_demeaning,
                                                                                Rcpp::Named("conv_demeaning") = conv_demeaning),
                                  Rcpp::Named("model_info") = Rcpp::List::create(Rcpp::Named("used_ids") = id_unique,
                                                                                 Rcpp::Named("beta_start") = beta_start,
                                                                                 Rcpp::Named("alpha_start") = alpha_start, 
                                                                                 Rcpp::Named("y") = y,
                                                                                 Rcpp::Named("X") = X,
                                                                                 Rcpp::Named("id") = id,
                                                                                 Rcpp::Named("time") = time,
                                                                                 Rcpp::Named("index_id") = index_id,
                                                                                 Rcpp::Named("Ti_vector") = Ti_vector));
      }
    case 1 :
      {
        // 1 - Analytical bias correction
        
        // Compute Uit, vit, and vit_alpha
        const arma::colvec mu = X * beta + alpha(index_id);
        const arma::colvec vit = gamma(y, mu, model);
        const arma::mat uit = X.each_col() % vit;
        const arma::colvec vit_alpha = - weights(y, mu, model);
        
        // Compute vit_sq, and V2it
        const arma::colvec vit_sq = arma::pow(vit, 2);
        const arma::colvec V2it = vit_sq + vit_alpha;
        
        // Compute sum of uit.each_col() % vit
        arma::mat sum_uv(nt, d);
        for (arma::uword j = 0 ; j < d ; ++j) {
          
          const arma::colvec uv = uit.col(j) % vit;
          start = 0;
          for (arma::uword i = 0 ; i < n ; ++i) {
            
            const arma::uword end = Ti_vector(i) + start - 1;
            
            sum_uv(arma::span(start, end), j).fill(arma::accu(uv.subvec(start, end)));
            
            start = end + 1;
          }
        }
        
        // Compute sum of vit_sq
        arma::colvec sum_vit_sq(nt);
        start = 0;
        for (arma::uword i = 0 ; i < n ; ++i) {
          
          const arma::uword end = Ti_vector(i) + start - 1;
          
          sum_vit_sq.subvec(start, end).fill(arma::accu(vit_sq.subvec(start, end)));
          
          start = end + 1;
        }
        
        // Compute Uit
        const arma::mat Uit = uit - vit % (sum_uv.each_col() / sum_vit_sq).each_col();
        
        // Compute B_bar, H_bar, and b_bar
        const arma::mat H_bar = Uit.t() * Uit / arma::mean(Ti_vector);
        const arma::colvec b_bar = Uit.t() * (V2it / (2.0 * sum_vit_sq));
        const arma::colvec B_bar = - H_bar.i() * b_bar;
        
        beta_tilde = beta - B_bar / arma::mean(Ti_vector);
        break;
      }
  }
  
  // Compute fixed effects after bias correction (Step 8 - pseudo-code)
  unsigned int iter_offset;
  bool conv_offset;
  const arma::colvec alpha_tilde = offset(iter_offset, conv_offset, beta_tilde, y, X, index_id, Ti_vector, model, iter_max2, tolerance2);
  
  // Compute standard errors after bias correction
  arma::colvec se_beta_tilde(d);
  arma::colvec se_alpha_tilde(n);
  arma::mat H_beta_tilde(d, d);
  standard_errors(se_beta_tilde, se_alpha_tilde, H_beta_tilde, beta_tilde, alpha_tilde, y, X, index_id, Ti_vector, model);
  
  
  // Return "1 - Analytical bias correction"
  const arma::colvec q = 2.0 * y - 1.0;
  return Rcpp::List::create(Rcpp::Named("par") = Rcpp::List::create(Rcpp::Named("beta") = beta,
                                                                    Rcpp::Named("alpha") = alpha,
                                                                    Rcpp::Named("se_beta") = se_beta,
                                                                    Rcpp::Named("se_alpha") = se_alpha,
                                                                    Rcpp::Named("beta_vcov") = H_beta,
                                                                    Rcpp::Named("avg_alpha") = arma::mean(alpha)),
                            Rcpp::Named("par_corr") = Rcpp::List::create(Rcpp::Named("beta") = beta_tilde,
                                                                         Rcpp::Named("alpha") = alpha_tilde,
                                                                         Rcpp::Named("se_beta") = se_beta_tilde,
                                                                         Rcpp::Named("se_alpha") = se_alpha_tilde,
                                                                         Rcpp::Named("beta_vcov") = H_beta_tilde,
                                                                         Rcpp::Named("avg_alpha") = arma::mean(alpha_tilde)),
                            Rcpp::Named("logl_info") = Rcpp::List::create(Rcpp::Named("loglik") = log_likelihood(q % (X * beta + alpha(index_id)), model),
                                                                          Rcpp::Named("iter_demeaning") = iter_demeaning,
                                                                          Rcpp::Named("conv_demeaning") = conv_demeaning,
                                                                          Rcpp::Named("loglik_corr") = log_likelihood(q % (X * beta_tilde + alpha_tilde(index_id)), model),
                                                                          Rcpp::Named("iter_offset") = iter_offset,
                                                                          Rcpp::Named("conv_offset") = conv_offset), 
                            Rcpp::Named("model_info") = Rcpp::List::create(Rcpp::Named("used_ids") = id_unique,
                                                                           Rcpp::Named("beta_start") = beta_start,
                                                                           Rcpp::Named("alpha_start") = alpha_start,
                                                                           Rcpp::Named("y") = y,
                                                                           Rcpp::Named("X") = X,
                                                                           Rcpp::Named("id") = id,
                                                                           Rcpp::Named("time") = time,
                                                                           Rcpp::Named("index_id") = index_id,
                                                                           Rcpp::Named("Ti_vector") = Ti_vector));
}
