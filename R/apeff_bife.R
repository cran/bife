#' @title
#' Average Partial Effects for Binary Choice Models with Fixed Effects
#' 
#' @description 
#' \code{apeff_bife} is a function used to compute average partial effects
#' for fixed effects binary choice models. It is able to compute
#' bias-corrected average partial effects derived by Newey and Hahn (2004) to
#' account for the incidental parameters bias.
#' 
#' @param 
#' mod an object of class \code{bife}.
#' 
#' @param 
#' discrete a description of the variables that are discrete regressors.
#' For \code{apeff_bife} this has to be a character string naming the discrete regressors.
#' Default is \code{NULL} (no discrete regressor(s)).
#' 
#' @param 
#' bias_corr an optional string that specifies the type of the bias correction:
#' semi or analytical. The value should be any of the values \code{"semi"} or 
#' \code{"ana"}. Default is \code{"ana"} (analytical bias-correction). Details
#' are given under \code{Details}.
#' 
#' @param 
#' iter_demeaning an optional integer value that specifies the maximum number
#' of iterations of the demeaning algorithm. Default is \code{100}. Details
#' are given under \code{Details}.
#' 
#' @param
#' tol_demeaning an optional number that specifies the tolerance level of the
#' demeaning algorithm. Default is \code{1e-5}. Details are given under \code{Details}.
#' 
#' @param 
#' iter_offset an optional integer value that specifies the maximum number of
#' iterations of the offset algorithm for the computation of bias-adjusted fixed
#' effects. Default is \code{1000}. Details are given under \code{Details}.
#' 
#' @param 
#' tol_offset an optional number that specifies the tolerance level of the 
#' offset algorithm for the computation of bias-adjusted fixed effects. Default
#' is \code{1e-5}. Details are given under \code{Details}.
#' 
#' @details
#' The semi bias-corrected average partial effects are computed as usual partial
#' effects with the bias-adjusted fixed effects and the bias-corrected structural
#' parameters.
#' 
#' The analytical bias-corrected average partial effects follow
#' Newey and Hahn (2004). For further details consult the description of \code{bife}.
#' 
#' \strong{Note:} Bias-corrected partial effects can be only returned if the
#' object \code{mod} returns bias-corrected coefficients, i.e. if a bias-correction 
#' has been used in the previous \code{bife} command.
#' 
#' @return 
#' An object of \code{apeff_bife} returns a named matrix with at least a first
#' column "apeff" containing the uncorrected average partial effects of the
#' structural variables. An optional second column "apeff_corrected" is returned
#' containing the corrected average partial effects of the structural variables.
#' 
#' @author
#' Amrei Stammann, Daniel Czarnowske, Florian Heiss, Daniel McFadden
#'
#' @seealso
#' \code{\link{bife}}
#' 
#' @references 
#' Hahn, J., and W. Newey (2004). "Jackknife and analytical bias reduction 
#' for nonlinear panel models." Econometrica 72(4), 1295-1319.
#' 
#' @references
#' Stammann, A., F. Heiss, and D. McFadden (2016). "Estimating Fixed Effects 
#' Logit Models with Large Panel Data". Working paper.
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
#' # Compute uncorrected average partial effects for mod_no
#' # Note: bias_corr does not affect the result
#' apeff_bife(mod_no, discrete = c("KID1", "KID2", "KID3")) 
#' 
#' # Fixed effects logit model with analytical bias-correction
#' mod_ana <- bife(LFP ~ AGE + I(INCH / 1000) + KID1 + KID2 + KID3 | ID,
#'  data = dataset)
#'                
#' # Compute semi-corrected average partial effects for mod_ana
#' apeff_bife(mod_ana, discrete = c("KID1", "KID2", "KID3"),
#'  bias_corr = "semi")
#' 
#' # Compute analytical bias-corrected average partial effects
#' # for mod_ana
#' apeff_bife(mod_ana, discrete = c("KID1", "KID2", "KID3"))
#'  
#' @export
apeff_bife <- function(mod, discrete = NULL,                       # object obtained from bife() and names of discrete variables 
                       bias_corr = "ana",                          # Type of bias correction
                       iter_demeaning = 100, tol_demeaning = 1e-5, # Control parameters "pseudo demeaning"
                       iter_offset = 1000, tol_offset = 1e-5) {    # Control parameters "glm offset"
 
 
  # Checking input arguments
  if(!inherits(mod, "bife")) stop("'apeff_bife' called on a non-'bife' object.")
  if(is.null(mod[["par_corr"]])) {
  
    corr <- "uncorrected"
  } else {
    
    if(bias_corr != "semi" && bias_corr != "ana") stop("'bias_corr' must be 'semi' or 'ana'.")
    corr <- "corrected"
  }
 
  # Name of structural parameters and map "discrete"
  X_name <- mod[["model_info"]][["str_name"]]
  discrete_vec <- numeric(length(X_name))
  discrete_vec[match(discrete, X_name)] = 1
 
  # Map "model"
  switch(mod[["model_info"]][["model"]], logit = model_int <- 0, probit = model_int <- 1)
 
  # Switch based on corr
  switch(corr,
         uncorrected = {
           apeff <- .avg_peff(discrete = discrete_vec, 
                              beta = mod[["par"]][["beta"]],
                              alpha = mod[["par"]][["alpha"]], 
                              X = mod[["model_info"]][["X"]],
                              index_id = mod[["model_info"]][["index_id"]],
                              model = model_int)
           colnames(apeff) <- "apeff"
           rownames(apeff) <- X_name},
         corrected = {
           switch(bias_corr, semi = bias_corr_int <- 0, ana = bias_corr_int <- 1)
           apeff <- .avg_peff_corr(discrete = discrete_vec, 
                                   beta = mod[["par"]][["beta"]],
                                   alpha = mod[["par"]][["alpha"]],
                                   beta_start = mod[["model_info"]][["beta_start"]],
                                   alpha_start = mod[["model_info"]][["alpha_start"]],
                                   beta_tilde = mod[["par_corr"]][["beta"]],
                                   alpha_tilde = mod[["par_corr"]][["alpha"]],
                                   y = mod[["model_info"]][["y"]],
                                   X = mod[["model_info"]][["X"]],
                                   time = mod[["model_info"]][["time"]],
                                   index_id = mod[["model_info"]][["index_id"]],
                                   Ti_vector = mod[["model_info"]][["Ti_vector"]], 
                                   model = model_int,
                                   bias_corr = bias_corr_int, 
                                   iter_max1 = iter_demeaning,
                                   tolerance1 = tol_demeaning,
                                   iter_max2 = iter_offset,
                                   tolerance2 = tol_offset)
           colnames(apeff) <- c("apeff", "apeff_corrected")
           rownames(apeff) <- X_name})
 
  # Adjustment due to perfect classification
  adjustment <- (mod[["logl_info"]][["nobs"]] - mod[["model_info"]][["drop_pc"]]) / mod[["logl_info"]][["nobs"]] 
 
 
  return(adjustment * apeff)
}
