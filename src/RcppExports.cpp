// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// avg_peff
arma::colvec avg_peff(const arma::colvec& discrete, const arma::colvec& beta, const arma::colvec& alpha, const arma::mat& X, const arma::uvec& index_id, const unsigned int model);
RcppExport SEXP _bife_avg_peff(SEXP discreteSEXP, SEXP betaSEXP, SEXP alphaSEXP, SEXP XSEXP, SEXP index_idSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type discrete(discreteSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type index_id(index_idSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(avg_peff(discrete, beta, alpha, X, index_id, model));
    return rcpp_result_gen;
END_RCPP
}
// avg_peff_corr
arma::mat avg_peff_corr(const arma::colvec& discrete, const arma::colvec& beta, const arma::colvec& alpha, const arma::colvec& beta_start, const arma::colvec& alpha_start, const arma::colvec& beta_tilde, const arma::colvec& alpha_tilde, const arma::colvec& y, const arma::mat& X, const arma::colvec& time, const arma::uvec& index_id, const arma::colvec& Ti_vector, const unsigned int model, const unsigned int bias_corr, const unsigned int iter_max1, const double tolerance1, const unsigned int iter_max2, const double tolerance2);
RcppExport SEXP _bife_avg_peff_corr(SEXP discreteSEXP, SEXP betaSEXP, SEXP alphaSEXP, SEXP beta_startSEXP, SEXP alpha_startSEXP, SEXP beta_tildeSEXP, SEXP alpha_tildeSEXP, SEXP ySEXP, SEXP XSEXP, SEXP timeSEXP, SEXP index_idSEXP, SEXP Ti_vectorSEXP, SEXP modelSEXP, SEXP bias_corrSEXP, SEXP iter_max1SEXP, SEXP tolerance1SEXP, SEXP iter_max2SEXP, SEXP tolerance2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type discrete(discreteSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta_start(beta_startSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type alpha_start(alpha_startSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta_tilde(beta_tildeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type alpha_tilde(alpha_tildeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type index_id(index_idSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Ti_vector(Ti_vectorSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type bias_corr(bias_corrSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type iter_max1(iter_max1SEXP);
    Rcpp::traits::input_parameter< const double >::type tolerance1(tolerance1SEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type iter_max2(iter_max2SEXP);
    Rcpp::traits::input_parameter< const double >::type tolerance2(tolerance2SEXP);
    rcpp_result_gen = Rcpp::wrap(avg_peff_corr(discrete, beta, alpha, beta_start, alpha_start, beta_tilde, alpha_tilde, y, X, time, index_id, Ti_vector, model, bias_corr, iter_max1, tolerance1, iter_max2, tolerance2));
    return rcpp_result_gen;
END_RCPP
}
// bife
Rcpp::List bife(const arma::colvec& y, const arma::mat& X, const arma::colvec& id, const arma::colvec& beta_start, const unsigned int model, const unsigned int bias_corr, const unsigned int iter_max1, const double tolerance1, const unsigned int iter_max2, const double tolerance2);
RcppExport SEXP _bife_bife(SEXP ySEXP, SEXP XSEXP, SEXP idSEXP, SEXP beta_startSEXP, SEXP modelSEXP, SEXP bias_corrSEXP, SEXP iter_max1SEXP, SEXP tolerance1SEXP, SEXP iter_max2SEXP, SEXP tolerance2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type id(idSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta_start(beta_startSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type bias_corr(bias_corrSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type iter_max1(iter_max1SEXP);
    Rcpp::traits::input_parameter< const double >::type tolerance1(tolerance1SEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type iter_max2(iter_max2SEXP);
    Rcpp::traits::input_parameter< const double >::type tolerance2(tolerance2SEXP);
    rcpp_result_gen = Rcpp::wrap(bife(y, X, id, beta_start, model, bias_corr, iter_max1, tolerance1, iter_max2, tolerance2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bife_avg_peff", (DL_FUNC) &_bife_avg_peff, 6},
    {"_bife_avg_peff_corr", (DL_FUNC) &_bife_avg_peff_corr, 18},
    {"_bife_bife", (DL_FUNC) &_bife_bife, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_bife(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
