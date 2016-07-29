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
 
 cat("Fixed effects", x$model.info$model, "model\n")
 
 switch(x$model.info$bias.corr,
        no = bias.corr <- "no",
        ana = bias.corr <- "analytical",
        jack = bias.corr <- "jackknife")
 
 cat("with", bias.corr, "bias-correction\n\n")
 
 cat("Uncorrected estimate(s)=", round(x$par$beta, digits = digits), "\n\n")
 
 if(!is.null(x$par.corr)) {
  
  cat("Corrected estimate(s)=", round(x$par.corr$beta, digits = digits), "\n\n")
 }
}