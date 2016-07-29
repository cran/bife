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
 if(is.null(object$par.corr)) corrected <- FALSE
 corr <- ifelse(corrected == TRUE, "corrected", "uncorrected")
 alpha <- ifelse(fixed == TRUE, "full", "average")
 
 df <- object$logl.info$df
 X.name <- object$model.info$str.name
 coef.colnames <- c("Estimate", "Std. error", "t-value", "Pr(> t)")
 
 
 switch(alpha,
        average = {switch(corr,
                          corrected = {coef.beta <- object$par.corr$beta
                                       se.beta <- object$par.corr$se.beta
                                       t.value.beta <- coef.beta / se.beta
                                       tabular.beta <- cbind(coef.beta, se.beta, t.value.beta, 2.0 * pt(-abs(t.value.beta), df = df))
                                       rownames(tabular.beta) <- X.name
                                       colnames(tabular.beta) <- coef.colnames
                                       results <- list(coefmat.beta = tabular.beta,
                                                       avg.alpha = object$par.corr$avg.alpha,
                                                       loglik = object$logl.info$loglik.corr,
                                                       iter.offset = object$logl.info$iter.offset,
                                                       conv.offset = object$logl.info$conv.offset)},
                          uncorrected = {coef.beta <- object$par$beta
                                         se.beta <- object$par$se.beta
                                         t.value.beta <- coef.beta / se.beta
                                         tabular.beta <- cbind(coef.beta, se.beta, t.value.beta, 2.0 * pt(-abs(t.value.beta), df = df))
                                         rownames(tabular.beta) <- X.name
                                         colnames(tabular.beta) <- coef.colnames
                                         results <- list(coefmat.beta = tabular.beta,
                                                         avg.alpha = object$par$avg.alpha,
                                                         loglik = object$logl.info$loglik)})},
        full = {switch(corr,
                       corrected = {coef.beta <- object$par.corr$beta
                                    se.beta <- object$par.corr$se.beta
                                    t.value.beta <- coef.beta / se.beta
                                    tabular.beta <- cbind(coef.beta, se.beta, t.value.beta, 2.0 * pt(-abs(t.value.beta), df = df))
                                    rownames(tabular.beta) <- X.name
                                    colnames(tabular.beta) <- coef.colnames
                                    coef.alpha <- object$par.corr$alpha
                                    se.alpha <- object$par.corr$se.alpha
                                    t.value.alpha <- coef.alpha / se.alpha
                                    tabular.alpha <- cbind(coef.alpha, se.alpha, t.value.alpha, 2.0 * pt(-abs(t.value.alpha), df = df))
                                    rownames(tabular.alpha) <- object$model.info$used.ids
                                    colnames(tabular.alpha) <- coef.colnames
                                    results <- list(coefmat.beta = tabular.beta,
                                                    coefmat.alpha = tabular.alpha,
                                                    loglik = object$logl.info$loglik.corr,
                                                    iter.offset = object$logl.info$iter.offset,
                                                    conv.offset = object$logl.info$conv.offset)},
                       uncorrected = {coef.beta <- object$par$beta
                                      se.beta <- object$par$se.beta
                                      t.value.beta <- coef.beta / se.beta
                                      tabular.beta <- cbind(coef.beta, se.beta, t.value.beta, 2.0 * pt(-abs(t.value.beta), df = object$logl.info$df))
                                      rownames(tabular.beta) <- X.name
                                      colnames(tabular.beta) <- coef.colnames
                                      coef.alpha <- object$par$alpha
                                      se.alpha <- object$par$se.alpha
                                      t.value.alpha <- coef.alpha / se.alpha
                                      tabular.alpha <- cbind(coef.alpha, se.alpha, t.value.alpha, 2.0 * pt(-abs(t.value.alpha), df = object$logl.info$df))
                                      rownames(tabular.alpha) <- object$model.info$used.ids
                                      colnames(tabular.alpha) <- coef.colnames
                                      results <- list(coefmat.beta = tabular.beta,
                                                      coefmat.alpha = tabular.alpha,
                                                      loglik = object$logl.info$loglik)})})
 
 results$nobs <- object$logl.info$nobs
 results$events <- object$logl.info$events
 results$iter.demeaning <- object$logl.info$iter.demeaning
 results$conv.demeaning <- object$logl.info$conv.demeaning
 results$formula <- object$model.info$formula
 results$model <- object$model.info$model
 results$drop.NA <- object$model.info$drop.NA
 results$drop.pc <- object$model.info$drop.pc
 results$bias.corr <- object$model.info$bias.corr
 results$corrected <- corrected
 
 
 return(structure(results, class = "summary.bife"))
}