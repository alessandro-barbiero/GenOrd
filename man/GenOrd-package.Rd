\name{GenOrd-package}
\alias{GenOrd-package}
\docType{package}
\title{
Simulation of Discrete Random Variables with Given Correlation Matrix and Marginal Distributions
}
\description{
The package implements a procedure for generating samples from a multivariate discrete random variable with pre-specified correlation matrix and marginal distributions. The marginal distributions are linked together through either a Gaussian or a Student's t copula.
The procedure is developed in two steps: the first step (function \code{\link{ordcont}}) sets up the Gaussian/t copula in order to achieve the desired correlation matrix on the target random discrete components; the second step (\code{\link{ordsample}}) generates samples from the target variables.
The procedure can handle both Pearson's and Spearman's correlations, and any finite support for the discrete variables.
The intermediate function \code{\link{contord}} computes the correlations of the multivariate discrete variable derived from correlated normal variables through discretization.
Function \code{\link{corrcheck}} returns the lower and upper bounds of the correlation coefficient of each pair of discrete variables given their marginal distributions, i.e., returns the range of feasible bivariate correlations.
Function \code{\link{estcontord}}, to be used with bivariate samples only, estimates the parameters of the Gaussian/t copula and possibly of the margins.

Compared to version 1.4.0, this version has introduced the multivariate Student's t copula as an alternative latent structure
}
\details{
\tabular{ll}{
Package: \tab GenOrd\cr
Type: \tab Package\cr
Version: \tab 2.0.0\cr
Date: \tab 2025-08-11\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
Alessandro Barbiero, Pier Alda Ferrari

Maintainer: Alessandro Barbiero <alessandro.barbiero@unimi.it>
}
\references{
P.A. Ferrari, A. Barbiero (2012) Simulating ordinal data, Multivariate Behavioral Research, 47(4), 566-589

A. Barbiero, P.A. Ferrari (2015) Simulation of correlated Poisson variables. Applied Stochastic Models in Business and Industry,  31(5), 669-680.

A. Barbiero, P.A. Ferrari (2017) An R package for the simulation of correlated discrete variables. Communications in Statistics-Simulation and Computation, 46(7), 5123-5140.
}
\keyword{ package }
\seealso{
\code{\link{contord}}, \code{\link{ordcont}}, \code{\link{corrcheck}}, \code{\link{ordsample}}, \code{\link{estcontord}}
}

