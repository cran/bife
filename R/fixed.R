#' @title
#' Additional fixed effects
#'
#' @description
#' \code{fixed} can be used to extend \code{bife} with additional fixed effects. 
#' The function generates the corresponding design matrix by expanding factors to a set of dummy variables. 
#' To avoid collinearity issues the last column of the design matrix is suppressed.
#'  
#' \strong{Note 1:} The bias-corrections offered by \code{bife} are not able to correct the incidential parameter bias of multiple fixed effects.
#' Thus the option \code{bias_corr = "no"} is mandatory. See Newey and Hahn (2004) for further details of the bias-corrections used by \code{bife}.
#' 
#' \strong{Note 2:} Since \code{fixed} generates a potentially large matrix, this approach might be limited by memory.
#'  
#' @param 
#' x the name of variable in the corresponding data frame. Can also be a vector.
#' 
#' @author
#' Amrei Stammann, Daniel Czarnowske, Florian Heiss, Daniel McFadden
#' 
#' @references 
#' Hahn, J., and W. Newey (2004). "Jackknife and analytical bias reduction for nonlinear panel models". Econometrica 72(4), 1295-1319.
#'
#' @examples
#' library("bife")
#' 
#' # Load 'psid' dataset
#' dataset <- psid
#' 
#' # Fixed effects logit model w/o bias-correction
#' # with additional 'Time' fixed effects
#' mod <- bife(LFP ~ AGE + I(INCH / 1000) + KID1 + KID2 + KID3 + fixed(TIME) | ID, 
#'  data = dataset, bias_corr = "no")
#'  
#' # Summary of structural parameters only
#' summary(mod)
#' 
#' @export
fixed <- function(x) {
  
  
  D <- model.matrix( ~ factor(x) + 0)
  D <- D[, - ncol(D)]
  str_names <- colnames(D)
  colnames(D) <- gsub(".*)", "_", str_names)
  
  
  return(D)
}
