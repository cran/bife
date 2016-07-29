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
 
 switch(x$bias.corr,
        no = bias.corr <- "no",
        ana = bias.corr <- "analytical",
        jack = bias.corr <- "jackknife"
 )
 
 cat("with", bias.corr, "bias-correction\n\n")
 
 cat("Estimated model:\n")
 print(x$formula)
 cat("\nLog-Likelihood=", round(x$loglik, digits = digits), "\n")
 cat("n= ", x$nobs, ", number of events= ", x$events, "\n", sep = "")
 cat("Demeaning", ifelse(x$conv.demeaning == TRUE, 
                         paste("converged after", x$iter.demeaning, "iteration(s)\n"), 
                         "did not converge\n"))
 if(x$bias.corr != "no") {
  
  cat("Offset", ifelse(x$conv.offset == TRUE,
                       paste("converged after", x$iter.offset, "iteration(s)\n"),
                       "did not converge\n"))
 }
 
 if(x$corrected == TRUE) {
  
  cat("\nCorrected structural parameter(s):\n\n")
 } else {
  
  cat("\nUncorrected structural parameter(s):\n\n")
 }
 printCoefmat(x$coefmat.beta, P.values = TRUE, has.Pvalue = TRUE, digits = digits)
 cat("\n")
 
 if(x$drop.NA > 0) cat("(", x$drop.NA, "observation(s) deleted due to missingness )\n")
 
 if(x$drop.pc > 0) cat("(", x$drop.pc, "observation(s) deleted due to perfect classification )\n")
 
 if(is.null(x$coefmat.alpha)) {
  
  cat("Average individual fixed effects=", round(x$avg.alpha, digits = digits))
 } else {
  
  cat("Indiviudal fixed effects:\n\n")
  printCoefmat(x$coefmat.alpha, P.values = TRUE, has.Pvalue = TRUE, digits = digits)
 }
 
 cat("\n---------------------------------------------------------------\n")
}