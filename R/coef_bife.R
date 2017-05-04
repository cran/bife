#' @title
#' Extract Model Coefficients
#' 
#' @description
#' \code{coef.bife} is a generic function which extracts
#' model coefficients from objects returned by \code{bife}.
#' 
#' @param 
#' object an object of class \code{bife}.
#'
#' @param 
#' corrected an optional logical flag that specifies whether bias-corrected or uncorrected
#' coefficients are displayed. Default is \code{TRUE} (bias-corrected).
#' 
#' @param 
#' fixed an optional logical flag that specifies whether the structural parameters or the fixed 
#' effects are displayed. Default is \code{FALSE} (structural parameters).
#' 
#' @param 
#' ... other arguments
#'
#' @return
#' The function \code{coef.bife} returns a named vector of coefficients.
#' 
#' @author
#' Amrei Stammann, Daniel Czarnowske, Florian Heiss, Daniel McFadden
#' 
#' @seealso
#' \code{\link{bife}}
#' 
#' @export
coef.bife <- function(object, corrected = TRUE, fixed = FALSE, ...) {
 
  
  if(!inherits(object, "bife")) stop("'coef.bife' called on a non-'bife' object")
  if(is.null(object[["par_corr"]])) corrected <- FALSE
  corr <- ifelse(corrected == TRUE, "corrected", "uncorrected")
 
  switch(corr,
         corrected = {
           if(fixed == FALSE) {
             
             result <- as.vector(object[["par_corr"]][["beta"]])
             names(result) <- object[["model_info"]][["str_name"]]
           } else {
             result <- as.vector(object[["par_corr"]][["alpha"]])
             names(result) <- object[["model_info"]][["used_ids"]]
           }},
         uncorrected = {
           if(fixed == FALSE) {
             
             result <- as.vector(object[["par"]][["beta"]])
             names(result) <- object[["model_info"]][["str_name"]]
           } else {
             
             result <- as.vector(object[["par"]][["alpha"]])
             names(result) <- object[["model_info"]][["used_ids"]]
           }})
 
 
 return(result)
}
