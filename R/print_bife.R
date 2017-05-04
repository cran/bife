#' @title
#' Print \code{bife}
#' 
#' @description
#' \code{print.bife} is a generic function which displays some minimal
#' information from objects returned by \code{bife}.
#' 
#' @param 
#' x an object of class \code{bife}.
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
#' @export
print.bife <- function(x, digits = max(3, getOption("digits") - 3), ...) {
 
  
  if(!inherits(x, "bife")) stop("'print.bife' called on a non-'bife' object")
 
  cat("Fixed effects", x[["model_info"]][["model"]], "model\n")
 
  switch(x[["model_info"]][["bias_corr"]], no = bias_corr <- "no", ana = bias_corr <- "analytical")
 
  cat("with", bias_corr, "bias-correction\n\n")
 
  cat("Uncorrected estimate(s)=", round(x[["par"]][["beta"]], digits = digits), "\n\n")
 
  if(!is.null(x[["par_corr"]])) cat("Corrected estimate(s)=", round(x[["par_corr"]][["beta"]], digits = digits), "\n\n")
}