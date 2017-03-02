#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP bife_avg_peff(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bife_avg_peff_corr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bife_bife(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"bife_avg_peff",      (DL_FUNC) &bife_avg_peff,       6},
  {"bife_avg_peff_corr", (DL_FUNC) &bife_avg_peff_corr, 18},
  {"bife_bife",          (DL_FUNC) &bife_bife,          10},
  {NULL, NULL, 0}
};

void R_init_bife(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
