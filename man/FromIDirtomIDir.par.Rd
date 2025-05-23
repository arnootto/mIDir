% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FromIDtoMID.par.R
\name{FromIDirtomIDir.par}
\alias{FromIDirtomIDir.par}
\title{Convert inverted Dirichlet parameters to mode-based inverted Dirichlet parameterization}
\usage{
FromIDirtomIDir.par(alpha)
}
\arguments{
\item{alpha}{Numeric vector of concentration parameters for the inverted Dirichlet distribution, all elements must be greater than 0.}
}
\value{
A list with two elements:
\describe{
  \item{theta}{Numeric vector of length \code{d-1}, the mode parameters.}
  \item{gamma}{Numeric scalar, the dispersion parameter.}
}
}
\description{
Converts the parameters of an inverted Dirichlet distribution to mode-based inverted Dirichlet parameters theta and gamma.
}
\details{
The inverted Dirichlet distribution with parameters \code{alpha} = (\eqn{\alpha_1, \alpha_2, \dots, \alpha_d}) is reparameterized into mode parameters \code{theta} = (\eqn{\theta_1, \theta_2, \dots, \theta_{d-1}}) and dispersion parameter \code{gamma} (\eqn{\gamma}).
}
\examples{
alpha <- c(3, 4, 5)
FromIDtoMID.par(alpha)

}
