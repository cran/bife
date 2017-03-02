#' @title
#' Print \code{summary.bife}
#' 
#' @description
#' \code{print.summary.bife} is a generic function which displays summary statistics
#' from objects returned by \code{summary.bife}.
#' 
#' @param 
#' x an object of class \code{summary.bife}.
#'
#' @param 
#' digits integer indicating the number of decimal places. Default is 
#' \code{max(3, getOption("digits") - 3)}.
#' 
#' @param 
#' ... other arguments
#' 
#' @author
#' Amrei Stammann, Daniel Czarnowske, Florian Heiss, Daniel McFadden
#' 
#' @seealso
#' \code{\link{bife}}
#' 
#' @importFrom 
#' stats printCoefmat
#' 
#' @export
print.summary.bife <- function(x, digits = max(3, getOption("digits") - 3), ...) {
 
 if(!inherits(x, "summary.bife")) stop("'print.summary.bife' called on a non-'summary.bife' object")
 
 cat("---------------------------------------------------------------\n")
 cat("Fixed effects", x$model, "model\n")
 
 switch(x$bias_corr,
        no = bias_corr <- "no",
        ana = bias_corr <- "analytical",
        jack = bias_corr <- "jackknife"
 )
 
 cat("with", bias_corr, "bias-correction\n\n")
 
 cat("Estimated model:\n")
 print(x$formula)
 cat("\nLog-Likelihood=", round(x$loglik, digits = digits), "\n")
 cat("n= ", x$nobs, ", number of events= ", x$events, "\n", sep = "")
 cat("Demeaning", ifelse(x$conv_demeaning == TRUE, 
                         paste("converged after", x$iter_demeaning, "iteration(s)\n"), 
                         "did not converge\n"))
 if(x$bias_corr != "no") {
  
  cat("Offset", ifelse(x$conv_offset == TRUE,
                       paste("converged after", x$iter_offset, "iteration(s)\n"),
                       "did not converge\n"))
 }
 
 if(x$corrected == TRUE) {
  
  cat("\nCorrected structural parameter(s):\n\n")
 } else {
  
  cat("\nUncorrected structural parameter(s):\n\n")
 }
 printCoefmat(x$coefmat_beta, P.values = TRUE, has.Pvalue = TRUE, digits = digits)
 cat("\n")
 
 if(x$drop_NA > 0) cat("(", x$drop_NA, "observation(s) deleted due to missingness )\n")
 
 if(x$drop_pc > 0) cat("(", x$drop_pc, "observation(s) deleted due to perfect classification )\n")
 
 if(is.null(x$coefmat_alpha)) {
  
  cat("Average individual fixed effects=", round(x$avg_alpha, digits = digits))
 } else {
  
  cat("Indiviudal fixed effects:\n\n")
  printCoefmat(x$coefmat_alpha, P.values = TRUE, has.Pvalue = TRUE, digits = digits)
 }
 
 cat("\n---------------------------------------------------------------\n")
}