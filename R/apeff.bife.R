#' @title
#' Average Partial Effects for Binary Choice Models with Fixed Effects
#' 
#' @description 
#' \code{apeff.bife} is a function used to compute average partial effects
#' for fixed effects binary choice models. It is able to compute
#' bias-corrected average partial effects derived by Newey and Hahn (2004) to
#' account for the incidental parameters bias.
#' 
#' @param 
#' mod an object of class \code{bife}.
#' 
#' @param 
#' discrete a description of the variables that are discrete regressors. For \code{apeff.bife} this has to be a character string naming the discrete
#' regressors.
#' 
#' @param 
#' bias.corr an optional string that specifies the type of the bias correction:
#' semi, analytical, jackknife. The value should be any of the values \code{"semi"}, 
#' \code{"ana"}, or \code{"jack"}. Default is \code{"semi"} (semi bias-correction). Details
#' are given under \code{Details}.
#' 
#' @param 
#' iter.demeaning an optional integer value that specifies the maximum number
#' of iterations of the demeaning algorithm. Default is \code{100}. Details
#' are given under \code{Details}.
#' 
#' @param
#' tol.demeaning an optional number that specifies the tolerance level of the
#' demeaning algorithm. Default is \code{1e-5}. Details are given under \code{Details}.
#' 
#' @param 
#' iter.offset an optional integer value that specifies the maximum number of
#' iterations of the offset algorithm for the computation of bias-adjusted fixed
#' effects. Default is \code{1000}. Details are given under \code{Details}.
#' 
#' @param 
#' tol.offset an optional number that specifies the tolerance level of the 
#' offset algorithm for the computation of bias-adjusted fixed effects. Default
#' is \code{1e-5}. Details are given under \code{Details}.
#' 
#' @details
#' The semi bias-corrected average partial effects are computed as usual partial
#' effects with the bias-adjusted fixed effects and the bias-corrected structural
#' parameters depending on which bias-correction (analytical or jackknife) has
#' been used in the previously conducted \code{bife} command, i.e. which kind of
#' bias-corrected coefficients are stored in the object \code{mod}.
#' 
#' The analytical and jackknife bias-corrected average partial effects follow
#' Newey and Hahn (2004). The jackknife bias-correction requires to repeatedly
#' calling \code{bife} (IWLS demeaning algorithm) and to estimate the fixed
#' effects model on a subset of the data where the t-th period is excluded.
#' For further details consult the description of \code{bife}.
#' 
#' \strong{Note:} Bias-corrected partial effects can be only returned if the
#' object \code{mod} returns bias-corrected coefficients, i.e. if a bias-correction
#' (analytical or jackknife) has been used in the previous \code{bife} command.
#' 
#' @return 
#' An object of \code{apeff.bife} returns a named matrix with at least a first
#' column "apeff" containing the uncorrected average partial effects of the
#' structural variables. An optional second column "apeff.corrected" is returned
#' containing the corrected average partial effects of the structural variables
#' given the choosen bias-correction (semi, analytical, or jackknife).
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
#' data.set <- psid
#' head(data.set)
#' 
#' # Fixed effects logit model w/o bias-correction
#' mod.no <- bife(LFP ~ AGE + I(INCH / 1000) + KID1 + KID2 + KID3 | ID,
#'  data = data.set, bias.corr = "no")
#'                
#' # Compute uncorrected average partial effects for mod.no
#' # Note: bias.corr does not affect the result
#' apeff.bife(mod.no, discrete = c("KID1", "KID2", "KID3")) 
#' 
#' # Fixed effects logit model with analytical bias-correction
#' mod.ana <- bife(LFP ~ AGE + I(INCH / 1000) + KID1 + KID2 + KID3 | ID,
#'  data = data.set)
#'                
#' # Compute semi-corrected average partial effects for mod.ana
#' apeff.bife(mod.ana, discrete = c("KID1", "KID2", "KID3"))
#' 
#' # Compute analytical bias-corrected average partial effects
#' # for mod.ana
#' apeff.bife(mod.ana, discrete = c("KID1", "KID2", "KID3"),
#'  bias.corr = "ana")


apeff.bife <- function(mod, discrete,                              # object obtained from bife() and names of discrete variables 
                       bias.corr = "semi",                         # Type of bias correction
                       iter.demeaning = 100, tol.demeaning = 1e-5, # Control parameters "pseudo demeaning"
                       iter.offset = 1000, tol.offset = 1e-5) {    # Control parameters "glm offset"
 
 
 # Checking input arguments
 if(!inherits(mod, "bife")) stop("'apeff.bife' called on a non-'bife' object")
 if(is.null(mod$par.corr)) {
  
  corr <- "uncorrected"
 } else {
  
  if(bias.corr != "semi" && bias.corr != "ana" && bias.corr != "jack") stop("'bias.corr' must be 'semi', 'ana', or 'jack'")
  corr <- "corrected"
 }
 
 # Name of structural parameters and map "discrete"
 X.name <- mod$model.info$str.name
 discrete.vec <- numeric(length(X.name))
 discrete.vec[match(discrete, X.name)] = 1
 
 # Map "model"
 switch(mod$model.info$model,
        logit = model.int <- 0,
        probit = model.int <- 1
 )
 
 # Switch based on corr
 switch(corr,
        uncorrected = {apeff <- .apeff(discrete = discrete.vec, 
                                       beta = mod$par$beta, alpha = mod$par$alpha, 
                                       X = mod$model.info$X, index_id = mod$model.info$index.id,
                                       model = model.int)
        colnames(apeff) <- "apeff"
        rownames(apeff) <- X.name},
        corrected = {switch(bias.corr,
                            semi = bias.corr.int <- 0,
                            ana = bias.corr.int <- 1,
                            jack = bias.corr.int <- 2
        )
         apeff <- .apeff.corr(discrete = discrete.vec, 
                              beta = mod$par$beta, alpha = mod$par$alpha, beta_start = mod$model.info$beta.start,
                              beta_tilde = mod$par.corr$beta, alpha_tilde = mod$par.corr$alpha,
                              y = mod$model.info$y, X = mod$model.info$X, t = mod$model.info$t, index_id = mod$model.info$index.id,
                              T_vector = mod$model.info$T.vector, 
                              model = model.int, bias_corr = bias.corr.int, 
                              iter_max1 = iter.demeaning, tolerance1 = tol.demeaning,
                              iter_max2 = iter.offset, tolerance2 = tol.offset)
         colnames(apeff) <- c("apeff", "apeff.corrected")
         rownames(apeff) <- X.name}
 )
 
 
 return(apeff)
}
