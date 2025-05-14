#' Contaminated Inverted Dirichlet Density
#'
#' Computes the density (or log-density) of the contaminated mode-parameterized inverted Dirichlet distribution.
#'
#' @param x Numeric matrix or data frame with \code{d} columns, where each row represents a vector of positive values. The number of columns must be one more than the length of \code{theta}.
#' @param theta Numeric vector of length \code{d-1}, the mode parameters (non-negative).
#' @param gamma Numeric scalar, the variance parameter for the good component (positive).
#' @param delta Numeric scalar, the mixing proportion (between 0 and 1).
#' @param eta Numeric scalar, the inflation parameter for the contaminated component (positive).
#' @param log Logical; if \code{TRUE}, returns the log-density. If \code{FALSE}, returns the density. Default is \code{FALSE}.
#'
#' @details The contaminated inverted Dirichlet distribution has density:
#' \deqn{
#' f(x; \theta, \gamma, \delta, \eta) = (1-\delta) \cdot f_{\text{mIDir}}(x; \theta, \gamma) + \delta \cdot f_{\text{MID}}(x; \theta, \eta \gamma)
#' }
#' where \eqn{f_{\text{mIDir}}(x; \theta, \gamma)} is the density from \code{dmIDir}.
#'
#' @return A numeric vector of density values (if \code{log = FALSE}) or log-density values (if \code{log = TRUE}), one for each row of \code{x}.
#'
#' @examples
#' x <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
#' theta <- c(0.2, 0.3)
#' gamma <- 0.5
#' delta <- 0.2
#' eta <- 2
#' dcmIDir(x, theta, gamma, delta, eta, log = FALSE)
#'
#' @export
#'
dcmIDir <- function(x, theta, gamma, delta, eta){
  dd <- delta*dmIDir(x = x, param=c(theta, gamma), log = FALSE) + (1-delta)*dmIDir(x = x, param = c(theta, eta*gamma), log = FALSE)
  return(dd)
}
