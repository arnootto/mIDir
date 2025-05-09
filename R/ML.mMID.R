#' Finite Mixture of Mode-Parameterized Inverted Dirichlet Distributions
#'
#' Performs maximum likelihood estimation for a fintie mixture of mode-parameterized inverted Dirichlet distributions using the Expectation-Maximization (EM) algorithm.
#'
#' @param X Numeric matrix with \code{d} columns and \code{n} rows, where each row represents a vector of positive values.
#' @param k Integer, the number of mixture components (must be positive).
#' @param initialization Character string specifying the initialization method for posterior probabilities. Options are \code{"mclust"}, \code{"kmeans"}, \code{"random.soft"} (random soft assignments), \code{"random.hard"} (random hard assignments), or \code{"manual"} (uses \code{start.z}). Default is \code{"kmeans"}.
#' @param seed Integer, seed for random number generation if \code{initialization} is \code{"random.soft"} or \code{"random.hard"}. Default is \code{NULL}.
#' @param start.z Numeric matrix of size \code{n x k}, initial posterior probabilities for \code{initialization = "manual"}. Can be soft or hard assignments. Default is \code{NULL}.
#' @param method Optimization method within the M-step. Default is \code{"BFGS"}. Other options include \code{"Nelder-Mead"}.
#' @param maxit Maximum number of iterations for \code{optim} within step. Default is 1000.
#' @param reltol Relative convergence tolerance for \code{optim} within step. Default is 1e-15.
#' @param threshold Numeric, convergence threshold for the log-likelihood difference. Default is 1e-4.
#' @param max.iter Maximum number of EM iterations. Default is 5000.
#' @param plotll Logical; if \code{TRUE}, plots the log-likelihood over iterations. Default is \code{TRUE}.
#'
#' @details The \code{ML.mMID} function fits a mixture of \code{k} mode-parameterized inverted Dirichlet distributions to the data \code{X} using the EM algorithm.
#'
#' @return An object of class \code{"mixIDir"} (list) containing:
#' \describe{
#'   \item{X}{The input matrix \code{X}.}
#'   \item{k}{Number of mixture components.}
#'   \item{n}{Number of observations.}
#'   \item{density}{Matrix of size \code{n x k}, the weighted density for each observation and component.}
#'   \item{z}{Matrix of size \code{n x k}, posterior probabilities for each observation and component.}
#'   \item{prior}{Vector of length \code{k}, estimated mixing proportions.}
#'   \item{theta}{Matrix of size \code{d x k}, estimated mode parameters for each component.}
#'   \item{gamma}{Vector of length \code{k}, estimated dispersion parameters for each component.}
#'   \item{group}{Vector of length \code{n}, cluster assignments based on maximum posterior probability.}
#'   \item{iter.stop}{Number of EM iterations performed.}
#'   \item{npar}{Number of parameters estimated: \code{(k-1) + k*d}.}
#'   \item{loglik}{Final log-likelihood.}
#'   \item{AIC}{Akaike Information Criterion.}
#'   \item{BIC}{Bayesian Information Criterion.}
#'   \item{call}{The matched call.}
#' }
#'
#' @import stats
#' @export
ML.mMID <- function(
    X,
    k,
    initialization="kmeans",
    seed = NULL,
    start.z = NULL,
    method = "BFGS",
    maxit = 1000,
    reltol = 1e-15,
    threshold = 10^-4,
    max.iter = 5000,
    plotll = TRUE
){

  if(any(is.na(X)))
    stop('No NAs allowed.')
  if(is.null(k))
    stop('k is NULL')

  n <- nrow(X) # number of units
  d <- ncol(X) # number of dimensions

  z <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))

  # Parameters

  prior <- array(NA,c(k),dimnames=list(paste("comp.",1:k,sep="")))
  theta <- array(NA,c(d,k),dimnames=list(1:(d),paste("comp.",1:k,sep="")))  # matrix(NA, nrow = (d-1), ncol = k)
  gamma <- array(NA,c(k),dimnames=list(paste("comp.",1:k,sep="")))

  # Distribution

  dens <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))

  # posteriors initialization --------------------------------------------------

  if(k>1){

    if(initialization=="mclust"){
      z <- Mclust(X, G = k, modelNames = "V")$z
    }

    if(initialization=="kmeans"){
      clusters  <- kmeans(x = X, centers = k)
      z         <- clas(clusters$cluster,k)
    }

    if(initialization=="random.soft"){
      if(!is.null(seed))
        set.seed(seed)
      z  <- array(runif(n*k),c(n,k)) # soft posterior probabilities (no-normalized) (n x k)
      z  <- z/rowSums(z)             # soft posterior probabilities (n x k)
    }

    if(initialization=="random.hard"){
      if(!is.null(seed))
        set.seed(seed)
      z  <- t(rmultinom(n, size = 1, prob = rep(1/k,k)))  # hard posterior probabilities (n x k)
    }

    if(initialization=="manual"){ # z.start can be both soft and hard initialization
      z  <- start.z      # posterior probabilities (n x k) no-normalized
    }

  }
  if(k==1)
    z <- matrix(rep(1,n), nrow = n, ncol = k)

  # initial parameters

  prior <- colMeans(z)
  for(j in 1:k){

    temp <- ML.mID(X = X,
                    weights = z[,j],
                    method = method,
                    theta.init = NULL,
                    gamma.init = NULL,
                    maxit = maxit,
                    reltol = reltol)

    theta[,j] <- temp$theta
    gamma[j]  <- temp$gamma

    # --------- #
    # densities #
    # --------- #

    densX    <- dMID(x = X,param = c(theta[,j], gamma[j]), log = FALSE)
    dens[,j] <- prior[j]*densX

  }

  # ------------------------------------- #
  # Global - Observed-data log-likelihood #
  # ------------------------------------- #

  iteration <- 1
  loglik    <- NULL

  llvalues          <- sum(log(rowSums(dens)))
  loglik[iteration] <- llvalues

  # EM algorithm --------------------------------------------------

  check <- 0
  while(check<1){

    cat("*")
    iteration <- iteration + 1

    # ++++++ #
    # E-Step #
    # ++++++ #

    z <- dens/matrix(rep(rowSums(dens),k),ncol=k)

    # ++++++ #
    # M-Step #
    # ++++++ #

    # priors #

    prior <- colMeans(z)

    # UD parameters #

    for(j in 1:k){

      temp <- ML.mID(X = X,
                      weights = z[,j],
                      method = method,
                      theta.init = theta[,j],
                      gamma.init = gamma[j],
                      maxit = maxit,
                      reltol = reltol)

      theta[,j] <- temp$theta
      gamma[j]  <- temp$gamma

      # --------- #
      # densities #
      # --------- #

      densX    <- dMID(x = X, param = c(theta[,j], gamma[j]), log = FALSE)
      dens[,j] <- prior[j]*densX
    }
    print(iteration)
    # ------------------------------------- #
    # Global - Observed-data log-likelihood #
    # ------------------------------------- #

    llvalues <- sum(log(rowSums(dens)))
    loglik[iteration] <- llvalues
    diff = loglik[iteration]-loglik[iteration-1]
    # ------------- #
    # Stopping Rule #
    # ------------- #

    # dif <- loglik[iteration] - loglik[iteration-1]
    # if(dif < 0)
    #   warning("The log-likelihood is not monotone")
    # if(dif < threshold)
    #   check <- 1
    if(iteration == max.iter | k == 1 | diff<threshold)
      check <- 1

  }

  finalloglik <- loglik[iteration]

  # ---------------------------------------------------------- #
  # plot to check the monotonic behavior of the log-likelihood #
  # ---------------------------------------------------------- #

  if(plotll){

    par(mai=c(0.84,0.8,0.012,0.004))
    par(las = 3)
    par(cex.axis=0.7)
    par(cex.lab=1.2)
    plot(1:iteration,loglik[1:iteration],type="l",axes = FALSE,xlab="iterations",ylab="log-likelihood",lwd=2)
    axis(1, at = 1:iteration,label = 1:iteration)
    axis(2)
    box(col = "black")

  }

  # The EM-algorithm is finished #

  # --------------------- #
  # Classification Matrix #
  # --------------------- #

  group <- apply(z,1,which.max)

  # -------------------- #
  # Number of parameters #
  # -------------------- #

  npar <- (k-1) + k*d

  # -------------------- #
  # Information criteria #
  # -------------------- #

  AIC <- 2*finalloglik - 2*npar
  BIC <- 2*finalloglik - npar*log(n)

  result <- list(
    X         = X,
    k         = k,
    n         = n,
    density   = dens,
    z         = z,
    prior     = prior,
    theta     = theta,
    gamma     = gamma,
    group     = group,
    iter.stop = iteration,
    npar      = npar,
    loglik    = finalloglik,
    AIC       = AIC,
    BIC       = BIC,
    call      = match.call()
  )
  class(result) <- "mixIDir"

  return(result)

}
