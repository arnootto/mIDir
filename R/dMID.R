#' Mode-Parameterized Inverted Dirichlet Density
#'
#' Computes the density (or log-density) of the mode-parameterized inverted Dirichlet distribution.
#'
#' @param x Numeric matrix or data frame with \code{d} columns, where each row represents a vector of positive values. The number of columns must be one more than the length of \code{param[1:p]}.
#' @param param Numeric vector of length \code{p+1}, where \code{p = d-1}. The first \code{p} elements are the mode parameters \code{theta} (positive), and the last element is the variance parameter \code{gamma} (positive).
#' @param log Logical; if \code{TRUE}, returns the log-density. If \code{FALSE}, returns the density. Default is \code{FALSE}.
#'
#' @return A numeric vector of densities values (if \code{log = FALSE}) or log-density values (if \code{log = TRUE}), one for each row of \code{x}.
#'
#' @examples
#' x <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
#' param <- c(0.2, 0.3, 0.5)  # theta1, theta2, gamma
#' dMID(x, param, log = FALSE)
#' dMID(x, param, log = TRUE)
#'
#' @export
dMID <- function(x,param,log=FALSE){
  x <- data.frame(x)
  ncol <- ncol(x)
  p <-length(param)-1
  d <- (lgamma(sum(1+(2+p+1/param[p+1])*param[1:p])+2+1/param[p+1]))-sum(lgamma(1+(2+p+1/param[p+1])*param[1:p]))-lgamma(2+1/param[p+1])+rowSums(as.matrix(log(x[,1:ncol,drop=F])) %*%diag(param[1:p]*(2+p+1/param[p+1])))-((2+p+1/param[p+1])+(2+p+1/param[p+1])*sum(param[1:p]))*log(1+rowSums(x[,1:ncol,drop=F]))
  if (log==FALSE){
    return(exp(d))}
  if (log==TRUE){
    return(d)
  }
}
