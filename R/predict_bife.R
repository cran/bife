#' @title
#' Computes Predicted Probabilities 
#' 
#' @description
#' Returns the predicted probabilities of an object returned by \code{bife}.
#' 
#' @param 
#' object an object of class \code{bife}.
#'
#' @param 
#' X_new a regressor matrix for predictions. If not supplied predictions are based on the
#' matrix returned by the object \code{bife}. See \code{Details}.
#' 
#' @param 
#' alpha_new a scalar or vector of fixed effects. If not supplied predictions are based on the
#' vector of fixed effects returned by \code{bife}. See \code{Details}.
#'
#' @param 
#' corrected an optional logical flag that specifies whether the predicted probabilities are
#' based on the bias-corrected/-adjusted parameters. Default is \code{TRUE} (bias-corrected).
#' 
#' @param 
#' ... other arguments
#' 
#' @details 
#' The regressor matrix returned by the object \code{bife} only includes individuals that 
#' were not dropped during fitting due to a non-varying response (perfect classification).
#' The predicted probabilities of those observations are equal to their response.
#' 
#' If \code{alpha_new} is supplied as a scalar each predicted probability is computed with
#' the same fixed effect. If \code{alpha_new} is supplied as a vector it has to be of same 
#' dimension as the corresponding regressor matrix.
#' 
#' @return
#' The function \code{predict.bife} returns a (named) vector of predicted probabilities.
#' 
#' @author
#' Amrei Stammann, Daniel Czarnowske, Florian Heiss, Daniel McFadden
#' 
#' @examples
#' library("bife")
#' 
#' # Load 'psid' dataset
#' dataset <- psid
#' head(dataset)
#' 
#' # Fixed effects logit model w/o bias-correction
#' mod_no <- bife(LFP ~ AGE + I(INCH / 1000) + KID1 + KID2 + KID3 | ID, 
#'  data = dataset, bias_corr = "no")
#'  
#' # Compute predicted probabilities based on the regressor matrix
#' # and fixed effects stored in 'mod_no'
#' prob <- predict(mod_no)
#' 
#' # Compute predicted probabilities based on the regressor matrix
#' # and all fixed effects set to zero
#' prob_zero <- predict(mod_no, alpha_new = 0.0)
#' 
#' @seealso
#' \code{\link{bife}}
#' 
#' @export
predict.bife <- function(object, X_new = NULL, alpha_new = NULL, corrected = TRUE, ...) {
  
  
  if (!inherits(object, "bife")) stop("'predict.bife' called on a non-'bife' object.")
  if (is.null(object[["par_corr"]])) corrected <- FALSE
  if (corrected == TRUE) corr <- "corrected" else corr <- "uncorrected"
  
  # Switch CDF
  switch(object[["model_info"]][["model"]], logit = prob <- plogis , probit = prob <- pnorm)
  
  # Check for out-of-sample prediction
  if (is.null(X_new)) X <- object[["model_info"]][["X"]] else X <- X_new
  if (is.null(alpha_new)) {
    
    alpha <- coef(object, corrected = corrected, fixed = TRUE)[object[["model_info"]][["index_id"]] + 1]
    if (length(alpha) != nrow(X)) stop("'alpha_new' of wrong dimension.")
  } else {
    
    if (length(alpha_new) == 1) {
      
      alpha <- rep(alpha_new, nrow(X))
    } else if (length(alpha_new) == nrow(X)) {
      
      alpha <- alpha_new
    } else stop("'alpha_new' of wrong dimension.")
  }
  
  # Compute predicted probabilities
  pred <- as.vector(prob(X %*% coef(object, corrected = corrected) + alpha))
  if (is.null(X_new) & is.null(alpha_new)) names(pred) <- object[["model_info"]][["used_ids"]][object[["model_info"]][["index_id"]] + 1]
  
  
  return(pred)
}