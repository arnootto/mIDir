#' Convert mode-parameterized inverted Dirichlet paramaeters to inverted Dirichlet Parameters
#'
#' Converts mode-based inverted Dirichlet parameters, theta and gamma, to the parameters of an inverted Dirichlet distribution.
#'
#' @param theta Numeric vector of mode parameters, all elements must be positive.
#' @param gamma Numeric scalar, the dispersion parameter, must be positive.
#'
#' @return A numeric vector of length \code{d}, the concentration parameters \code{alpha}.
#'
#' @examples
#' theta <- c(0.2, 0.3)
#' gamma <- 0.5
#' FromMIDtoID.par(theta, gamma)
#'
#' @export
FromMIDtoID.par <- function(theta,gamma){

d <- length(theta) + 1
alpha <- numeric(d)
theta.plus <- sum(theta)
for(j in 1:(d-1))
  alpha[j] <- (1 +(2+(d-1)+1/gamma)*theta[j])
alpha[d] <- 2+1/gamma

return(alpha)

}
