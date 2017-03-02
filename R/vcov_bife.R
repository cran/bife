#' @title
#' Extract Covariance Matrix of Structural Parameters
#' 
#' @description
#' \code{vcov.bife} extracts the covariance matrix of the
#' structural paramters from objects returned by \code{bife}
#' 
#' @param 
#' object an object of class \code{bife}.
#'
#' @param 
#' corrected an optional logical flag that specifies whether the covariance
#' matrix of the bias-corrected or uncorrected structural parameters
#' are displayed. Default is \code{TRUE} (bias-corrected).
#' 
#' @return
#' The function \code{vcov.bife} returns a named covariance matrix of the
#' structural parameters.
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
vcov.bife <- function(object, corrected = TRUE, ...) {
 
 if(!inherits(object, "bife")) stop("'coef.bife' called on a non-'bife' object")
 if(is.null(object$par.corr)) corrected <- FALSE
 corr <- ifelse(corrected == TRUE, "corrected", "uncorrected")
 
 switch (corr,
         corrected = {result <- object$par_corr$beta_vcov
                      rownames(result) <- colnames(result) <- object$model_info$str_name},
         uncorrected = {result <- object$par$beta_vcov
                        rownames(result) <- colnames(result) <- object$model_info$str_name})
 
 
 return(result)
}