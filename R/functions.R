
###################
## Observed Mode ##
###################

moda <- function(X, weights = NULL){

  n <- nrow(X)
  d <- ncol(X)

  if(is.null(weights))
    weights = rep(1,n)/n
  if(sum(weights) != 1)
    weights = weights/sum(weights)

  theta <- numeric(d)
  for(j in 1:(d)){
    dd <- density(X[,j], weights = weights)
    theta[j] <- dd$x[which.max(dd$y)]
  }

  return(theta)

}

###################
## Initial gamma ##
###################

library(modi)
init.gamma <- function(X, weights = NULL, j = 1, theta){

  Xj <- X[,j]
  thetaj <- theta[j]
  n <- length(Xj)
  d <- ncol(X)

  if(is.null(weights))
    weights = rep(1,n)

  sigma2 <- modi::weighted.var(x = Xj, w = weights, na.rm = FALSE)
  f <- function(gamma, sigma2, thetaj, d){

    num <- gamma*(thetaj + gamma)*(1 - thetaj + d*gamma)
    den <- (1 + d*gamma)^2*(1 + (d + 1)*gamma)

    return(num/den - sigma2)

  }

  res <- uniroot(f = f, interval = c(0, 100), sigma2 = sigma2, thetaj = thetaj, d = d)

  return(res$root)

}
