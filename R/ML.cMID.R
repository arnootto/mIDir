#' Maximum Likelihood Estimation for Contaminated Inverted Dirichlet
#'
#' Performs maximum likelihood estimation for the contaminated inverted Dirichlet distribution.
#'
#' @param X Numeric matrix with \code{d} columns and \code{n} rows, where each row represents a vector of positive values.
#' @param weights Numeric vector of length \code{n}, weights for each observation. Defaults to a vector of 1's.
#' @param method Optimization method for \code{optim}. Default is "Nelder-Mead". Other options include "BFGS".
#' @param theta.init Numeric vector of length \code{d}, initial values for the mode parameters \code{theta}. If \code{NULL}, initialized using the empirical mode of X.
#' @param gamma.init Numeric scalar, initial value for the dispersion parameter \code{gamma}. If \code{NULL}, initialized using the emprical variance based on a method of moments function.
#' @param delta.init Numeric scalar, initial value for the mixing proportion \code{delta} (between 0 and 1). Defaults to 0.01 if \code{NULL}.
#' @param eta.init Numeric scalar, initial value for the inflation parameter \code{eta} (greater than 1). Defaults to 1.01 if \code{NULL}.
#' @param maxit Maximum number of iterations for \code{optim}. Default is 10000.
#' @param reltol Relative convergence tolerance for \code{optim}. Default is 1e-15.
#' @param trace Non-negative integer; if positive, tracing information is printed during optimization. Default is 0.
#'
#'
#' @return A list containing:
#' \describe{
#'   \item{X}{The input matrix \code{X}.}
#'   \item{good.prob}{Numeric vector of length \code{n}, the probability that each observation belongs to the good component.}
#'   \item{good}{Binary vector of length \code{n}, indicating whether each observation is classified as good (1) or contaminated (0) based on \code{good.prob}.}
#'   \item{theta}{Estimated mode parameters (vector of length \code{d}).}
#'   \item{gamma}{Estimated dispersion parameter (scalar).}
#'   \item{delta}{Estimated mixing proportion (scalar, between 0 and 1).}
#'   \item{eta}{Estimated inflation parameter (scalar, greater than 1).}
#'   \item{npar}{Number of parameters estimated (\code{d + 3}).}
#'   \item{loglik}{The log-likelihood at the maximum.}
#'   \item{AIC}{Akaike Information Criterion}.
#'   \item{BIC}{Bayesian Information Criterion}.
#'   \item{convergence}{Convergence code from \code{optim} (0 indicates successful convergence).}
#'   \item{counts}{Number of function and gradient evaluations from \code{optim}.}
#' }
#'
#' @examples
#' set.seed(123)
#' X <- matrix(exp(rnorm(30)), nrow = 10, ncol = 3)
#' wML.cMID(X, weights = rep(1, 10), theta.init = c(0.2, 0.3, 0.4),
#'          gamma.init = 0.5, delta.init = 0.8, eta.init = 2)
#'
#' @import stats
#' @export
ML.cMID <- function(X,
                     weights = NULL,
                     method = "Nelder-Mead",
                     theta.init = NULL,
                     gamma.init = NULL,
                     delta.init = NULL,
                     eta.init = NULL,
                     maxit = 10000,
                     reltol = 1e-15,
                     trace = 0){

  X  <- as.matrix(X)
  d  <- ncol(X)
  n  <- nrow(X)

  if(is.null(weights))
    weights <- rep(1,n)

  npar <- d + 3

  # -------------- #
  # Initialization #
  # -------------- #

  # true world

  if(is.null(theta.init)){
    # w.mat <- matrix(rep(weights,d),nrow=n,ncol=d,byrow=FALSE)
    # theta <- colMeans(w.mat*X)[1:d]
    theta.init <- moda(X = X, weights = weights)
  }
  theta <- theta.init
  if(is.null(gamma.init))
    gamma.init <- init.gamma(X = X, weights = weights, j = 1, theta = theta)
  gamma <- gamma.init
  #print(gamma)
  if(is.null(delta.init))
    delta.init <- 0.99
  delta  <- delta.init
  if(is.null(eta.init))
    eta.init <- 1.01
  eta   <- eta.init

  # optim world

  init <- numeric(npar)
  for(j in 1:(d))
    init[j] <- log(theta[j])
  init[d+1]   <- log(gamma)
  init[d+2] <- log(delta/(1-delta))
  init[d+3] <- log(eta - 1)

  # objective function

  f <- function(par, X, weights, d){

    theta <- numeric(d)
    ub <- 0
    for(j in 1:(d)){
      #prob <- exp(par[j])/(1+exp(par[j]))
      #theta[j] <- prob/(1-ub)
      theta[j] <- exp(par[j])
      ub <- ub + theta[j]
    }
    gamma <- exp(par[d+1])
    delta  <- (0.5 + exp(par[d+2]))/(1 + exp(par[d+2]))
    eta   <- exp(par[d+3]) + 1

    # -------------- #
    # Log-Likelihood #
    # -------------- #

    loglik <- sum(weights*log(dcMID(x = X, theta = theta, gamma = gamma, delta = delta, eta = eta)))

    return(loglik)

  }

  res <- optim(par = init, fn = f, X = X, weights = weights, d = d,
               method = method, control = list(fnscale = -1, maxit = maxit, reltol = reltol, trace = trace))

  loglik <- res$value
  par   <- res$par

  # estimated parameters

  theta <- numeric(d)
  ub <- 0
  for(j in 1:(d)){
    prob <- exp(par[j])#/(1+exp(par[j]))
    theta[j] <- prob
    ub <- ub + theta[j]
  }
  gamma <- exp(par[d+1])
  delta  <- (0.5 + exp(par[d+2]))/(1 + exp(par[d+2]))
  eta   <- exp(par[d+3]) + 1

  # ------------------------------------ #
  # Probability to be a good observation #
  # ------------------------------------ #

  good.prob <- delta*dMID(x = X, param = c(theta , gamma), log = FALSE)/dcMID(x = X, theta = theta, gamma = gamma, delta = delta, eta = eta)
  good <- round(good.prob,0)

  # --------------- #
  # Model Selection #
  # --------------- #

  AIC    <- 2*loglik - 2*npar
  BIC    <- 2*loglik - npar*log(n)

  return(
    list(
      X      = X,
      good.prob = good.prob,
      good   = good,
      theta  = theta,
      gamma  = gamma,
      delta   = delta,
      eta    = eta,
      npar   = npar,
      loglik = loglik,
      AIC    = AIC,
      BIC    = BIC,
      convergence = res$convergence, # see optim()
      counts      = res$counts # see optim()
    )
  )

}
