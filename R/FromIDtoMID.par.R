#' Convert inverted Dirichlet parameters to mode-based inverted Dirichlet parameterization
#'
#' Converts the parameters of an inverted Dirichlet distribution to mode-based inverted Dirichlet parameters theta and gamma.
#'
#' @param alpha Numeric vector of concentration parameters for the inverted Dirichlet distribution, all elements must be greater than 0.
#'
#' @details The inverted Dirichlet distribution with parameters \code{alpha} = (\eqn{\alpha_1, \alpha_2, \dots, \alpha_d}) is reparameterized into mode parameters \code{theta} = (\eqn{\theta_1, \theta_2, \dots, \theta_{d-1}}) and dispersion parameter \code{gamma} (\eqn{\gamma}).
##'
#' @return A list with two elements:
#' \describe{
#'   \item{theta}{Numeric vector of length \code{d-1}, the mode parameters.}
#'   \item{gamma}{Numeric scalar, the dispersion parameter.}
#' }
#'
#' @examples
#' alpha <- c(3, 4, 5)
#' FromIDtoMID.par(alpha)
#'
#' @export
FromIDtoMID.par <- function(alpha){
  d <- length(alpha)
  theta <- numeric(d-1)
  for(j in 1:(d-1))
    theta[j] <- (alpha[j]-1)/(alpha[d]+(d-1))
  gamma <- 1/(alpha[d]-2)
  return(
    list(
      theta = theta,
      gamma = gamma
    )
  )

}
