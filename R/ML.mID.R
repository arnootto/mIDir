#' Maximum Likelihood Estimation for Mode-Parameterized Inverted Dirichlet distribution
#'
#' Performs maximum likelihood estimation for the mode-parameterized inverted Dirichlet distribution.
#'
#' @param X Numeric matrix with \code{d} columns and \code{n} rows, where each row represents a vector of positive values.
#' @param weights Numeric vector of length \code{n}, weights for each observation. Defaults to a vector of 1's.
#' @param method Optimization method for \code{optim}. Default is "BFGS". Other options include "Nelder-Mead".
#' @param theta.init Numeric vector of length \code{d-1}, initial values for the mode parameters \code{theta}. If \code{NULL}, initialized as empirical mode.
#' @param gamma.init Numeric scalar, initial value for the dispersion parameter \code{gamma}. If \code{NULL}, defaults to 0.5.
#' @param maxit Maximum number of iterations for \code{optim}. Default is 1000.
#' @param reltol Relative convergence tolerance for \code{optim}. Default is 1e-15.
#'
#'
#' @return A list containing:
#' \describe{
#'   \item{X}{The input matrix \code{X}.}
#'   \item{theta}{Estimated mode parameters (vector of length \code{d-1}).}
#'   \item{gamma}{Estimated dispersion parameter (scalar).}
#'   \item{npar}{Number of parameters estimated (\code{d}).}
#'   \item{loglik}{The log-likelihood at the maximum.}
#'   \item{AIC}{Akaike Information Criterion.}
#'   \item{BIC}{Bayesian Information Criterion.}
#'   \item{convergence}{Convergence code from \code{optim} (0 indicates successful convergence).}
#'   \item{counts}{Number of function and gradient evaluations from \code{optim}.}
#' }
#'
#' @import stats
#' @export
ml.mIDir <- function(X,
                    weights = NULL,
                    method = "BFGS",
                    theta.init = NULL,
                    gamma.init = NULL,
                    maxit = 1000,
                    reltol = 1e-15){

  X  <- as.matrix(X)
  d  <- ncol(X)
  n  <- nrow(X)

  if(is.null(weights))
    weights <- rep(1,n)

  # -------------- #
  # Initialization #
  # -------------- #

  if(is.null(theta.init)){

    theta.init <- moda(X = X, weights = weights)
  }
  if(is.null(gamma.init))
    gamma.init <- 0.5

  theta <- log(theta.init)
  gamma <- log(gamma.init)


  # initialization

  init <- c(theta,gamma)

  # objective function

  f <- function(par, X, weights){
    # p=length(par)-1
    # theta <- exp(par[1:p])
    # gamma <- exp(par[p+1])
    par=exp(par)
    # -------------- #
    # Log-Likelihood #
    # -------------- #

    loglik <- sum(weights*dmIDir(x = X, param=par, log = TRUE))

    return(loglik)

  }

  res <- optim(par = init, fn = f, X = X, weights = weights,
               method = method, control = list(fnscale = -1, maxit = maxit, reltol = reltol))

  loglik <- res$value
  par    <- res$par

  p=length(par)-1
  theta <- exp(par[1:p])
  gamma <- exp(par[p+1])
  npar   <- d+1
  AIC    <- 2*loglik - 2*npar
  BIC    <- 2*loglik - npar*log(n)

  return(
    list(
      X      = X,
      theta  = theta,
      gamma  = gamma,
      #alpha  = alpha,
      npar   = npar,
      loglik = loglik,
      AIC    = AIC,
      BIC    = BIC)
  )

}

