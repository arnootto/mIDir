#' Kernel Density Estimation for Inverted Dirichlet Distribution
#'
#' Computes a kernel density estimate for a sample using the mode-parameterized inverted Dirichlet distribution as the kernel.
#'
#' @param x Numeric vector or matrix with \code{d} columns, the evaluation point(s) where the density is computed. If a vector, it is converted to a single-row matrix.
#' @param data Numeric matrix with \code{d} columns and \code{n} rows, the sample data where each row represents a vector of positive values.
#' @param gamma Numeric scalar, the variance parameter for the inverted Dirichlet kernel (positive).
#'
#' @details The \code{kernelID} function computes a kernel density estimate at point \code{x} using the mode-parameterized inverted Dirichlet distribution as the kernel.
#'
#' @return A numeric scalar, the kernel density estimate at \code{x}.
#'
#' @export
kernelID <- function(x,data,gamma){
  x <- matrix(x, nrow = 1)

  n <- nrow(data)
  d <- ncol(data)
  bump <- numeric(n)
  for(i in 1:n){
    bump[i] <- dMID(x = x, param = c(data[i,], gamma), log = FALSE)
    print(i)
  }

  return(mean(bump))
}
