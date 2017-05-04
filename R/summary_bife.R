#' @title
#' Summarizing Binary Choice Models with Fixed Effects
#' 
#' @description
#' Summary statistics for objects of class \code{bife}.
#'
#' @param 
#' object an object of class \code{bife}.
#'
#' @param
#' corrected an optional logical flag that specifies 
#' whether bias-corrected or uncorrected coefficients 
#' are displayed. Default is \code{TRUE} (bias-corrected).
#'
#' @param
#' fixed an optional logical flag that specifies whether 
#' only structural parameters or all coefficients
#' (structural parameters and fixed effects) are displayed.
#' Default is \code{FALSE} (only structural parameters).
#' 
#' @param 
#' ... other arguments
#'
#' @return
#' Returns an object of class \code{summary.bife} which is a list of summary
#' statistics of \code{object}.
#' 
#' @author
#' Amrei Stammann, Daniel Czarnowske, Florian Heiss, Daniel McFadden
#' 
#' @seealso
#' \code{\link{bife}}
#' 
#' @importFrom
#' stats pt
#' 
#' @export
summary.bife <- function(object, corrected = TRUE, fixed = FALSE, ...) {
 
  
  if(!inherits(object, "bife")) stop("'summary.bife' called on a non-'bife' object")
  if(is.null(object[["par_corr"]])) corrected <- FALSE
  corr <- ifelse(corrected == TRUE, "corrected", "uncorrected")
  alpha <- ifelse(fixed == TRUE, "full", "average")
 
  # Compute coefmat
  nt <- object[["logl_info"]][["nobs"]]
  k <- object[["logl_info"]][["k"]] 
  df <- nt - k
  X_name <- object[["model_info"]][["str_name"]]
  coef_colnames <- c("Estimate", "Std. error", "t-value", "Pr(> t)")
  switch(alpha,
         average = {switch(corr,
                           corrected = {coef_beta <- object[["par_corr"]][["beta"]]
                                        se_beta <- object[["par_corr"]][["se_beta"]]
                                        t_value_beta <- coef_beta / se_beta
                                        tabular_beta <- cbind(coef_beta, se_beta, t_value_beta, 2.0 * pt(- abs(t_value_beta), df = df))
                                        rownames(tabular_beta) <- X_name
                                        colnames(tabular_beta) <- coef_colnames
                                        results <- list(coefmat_beta = tabular_beta,
                                                        avg_alpha = object[["par_corr"]][["avg_alpha"]],
                                                        loglik = object[["logl_info"]][["loglik_corr"]],
                                                        iter_offset = object[["logl_info"]][["iter_offset"]],
                                                        conv_offset = object[["logl_info"]][["conv_offset"]])},
                           uncorrected = {coef_beta <- object[["par"]][["beta"]]
                                          se_beta <- object[["par"]][["se_beta"]]
                                          t_value_beta <- coef_beta / se_beta
                                          tabular_beta <- cbind(coef_beta, se_beta, t_value_beta, 2.0 * pt(- abs(t_value_beta), df = df))
                                          rownames(tabular_beta) <- X_name
                                          colnames(tabular_beta) <- coef_colnames
                                          results <- list(coefmat_beta = tabular_beta,
                                                          avg_alpha = object[["par"]][["avg_alpha"]],
                                                          loglik = object[["logl_info"]][["loglik"]])})},
         full = {switch(corr,
                        corrected = {coef_beta <- object[["par_corr"]][["beta"]]
                                     se_beta <- object[["par_corr"]][["se_beta"]]
                                     t_value_beta <- coef_beta / se_beta
                                     tabular_beta <- cbind(coef_beta, se_beta, t_value_beta, 2.0 * pt(- abs(t_value_beta), df = df))
                                     rownames(tabular_beta) <- X_name
                                     colnames(tabular_beta) <- coef_colnames
                                     coef_alpha <- object[["par_corr"]][["alpha"]]
                                     se_alpha <- object[["par_corr"]][["se_alpha"]]
                                     t_value_alpha <- coef_alpha / se_alpha
                                     tabular_alpha <- cbind(coef_alpha, se_alpha, t_value_alpha, 2.0 * pt(- abs(t_value_alpha), df = df))
                                     rownames(tabular_alpha) <- object[["model_info"]][["used_ids"]]
                                     colnames(tabular_alpha) <- coef_colnames
                                     results <- list(coefmat_beta = tabular_beta,
                                                     coefmat_alpha = tabular_alpha,
                                                     loglik = object[["logl_info"]][["loglik_corr"]],
                                                     iter_offset = object[["logl_info"]][["iter_offset"]],
                                                     conv_offset = object[["logl_info"]][["conv_offset"]])},
                        uncorrected = {coef_beta <- object[["par"]][["beta"]]
                                       se_beta <- object[["par"]][["se_beta"]]
                                       t_value_beta <- coef_beta / se_beta
                                       tabular_beta <- cbind(coef_beta, se_beta, t_value_beta, 2.0 * pt(- abs(t_value_beta), df = df))
                                       rownames(tabular_beta) <- X_name
                                       colnames(tabular_beta) <- coef_colnames
                                       coef_alpha <- object[["par"]][["alpha"]]
                                       se_alpha <- object[["par"]][["se_alpha"]]
                                       t_value_alpha <- coef_alpha / se_alpha
                                       tabular_alpha <- cbind(coef_alpha, se_alpha, t_value_alpha, 2.0 * pt(- abs(t_value_alpha), df = df))
                                       rownames(tabular_alpha) <- object[["model_info"]][["used_ids"]]
                                       colnames(tabular_alpha) <- coef_colnames
                                       results <- list(coefmat_beta = tabular_beta,
                                                       coefmat_alpha = tabular_alpha,
                                                       loglik = object[["logl_info"]][["loglik"]])})})
 
  # Add information to results
  results[["nobs"]] <- object[["logl_info"]][["nobs"]]
  results[["events"]] <- object[["logl_info"]][["events"]]
  results[["iter_demeaning"]] <- object[["logl_info"]][["iter_demeaning"]]
  results[["conv_demeaning"]] <- object[["logl_info"]][["conv_demeaning"]]
  results[["formula"]] <- object[["model_info"]][["formula"]]
  results[["model"]] <- object[["model_info"]][["model"]]
  results[["drop_NA"]] <- object[["model_info"]][["drop_NA"]]
  results[["drop_pc"]] <- object[["model_info"]][["drop_pc"]]
  results[["bias_corr"]] <- object[["model_info"]][["bias_corr"]]
  results[["corrected"]] <- corrected
  
  # Compute AIC and BIC
  results[["AIC"]] <- - 2.0 * results[["loglik"]] + 2.0 * k
  results[["BIC"]] <- - 2.0 * results[["loglik"]] + k * log(nt)
  
  
  return(structure(results, class = "summary.bife"))
}