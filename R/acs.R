#' @title 
#' Female labor force participation - "American Community Survey (ACS PUMS 2014)"
#' 
#' @description 
#' The sample is drawn from the American Community Survey (ACS PUMS 2014) were the
#' panel structure is slightly different in comparison to the "classic" structure.
#' Overall 662,775 married women in \eqn{N = 51} states were observed. Since
#' each state is of different population size, this results in a highly
#' unbalanced panel were the largest state consists of \eqn{T_{max} = 74,752}
#' and the smallest of \eqn{T_{min} = 855} married women.
#'
#' @format A data frame with 662,775 rows:
#' \describe{
#'   \item{ST}{state identifier}
#'   \item{AGEP}{age of woman}
#'   \item{FER}{indicates if a woman gave birth to a child within the past 12 months}
#'   \item{PINCP}{total persons income}
#'   \item{LFP}{labor force participation}}
#'   
#' @references
#' American Community Survey. \url{http://www.census.gov}.
#'   
#' @seealso
#' \code{\link{bife}}
"acs"