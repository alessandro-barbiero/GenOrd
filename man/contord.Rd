\name{contord}
\alias{contord}
\title{
Correlations of discretized variables
}
\description{
The function computes the correlation matrix of the \eqn{k} variables, with given marginal distributions, derived discretizing a \eqn{k}-variate standard normal or Student's \eqn{t} variable with given correlation matrix
}
\usage{
contord(marginal, Sigma, support = list(), Spearman = FALSE, df=Inf, integerdf=TRUE, prob=FALSE)
}
\arguments{
  \item{marginal}{
a list of \eqn{k} elements, where \eqn{k} is the number of variables.
The \eqn{i}-th element of \code{marginal} is the vector of the cumulative probabilities defining the marginal distribution of the \eqn{i}-th component of the  multivariate variable. If the \eqn{i}-th component can take \eqn{k_i} values, the \eqn{i}-th element of \code{marginal} will contain \eqn{k_i-1} probabilities (the \eqn{k_i}-th is obviously 1 and shall not be included).
}
  \item{Sigma}{
the correlation matrix of the standard multivariate normal variable
}
  \item{support}{
a list of \eqn{k} elements, where \eqn{k} is the number of variables. The \eqn{i}-th element of \code{support} is the vector containing the ordered values of the support of the \eqn{i}-th variable. By default, the support of the \eqn{i}-th variable is \eqn{1,2,...,k_i}
}
  \item{Spearman}{
if \code{TRUE}, the function finds Spearman's correlations (and it is not necessary to provide \code{support}), if \code{FALSE} (default) Pearson's correlations
}
  \item{df}{
the degrees of freedom of the multivariate Student's \eqn{t}
}
  \item{integerdf}{
if \code{TRUE} uses \code{pmvt}, which requires an integer value for \code{df}; otherwise, uses \code{pmvt.alt}, which integrates \code{dmvt} with possibly non-integer \code{df}
}
  \item{prob}{
if \code{TRUE}, it return also the probability table of the bivariate discrete random variable
}
}
\value{
the correlation matrix of the discretized variables; if \code{prob=TRUE}, it returns a list containing the correlation along with the probability table, in case of a bivariate random variable
}
\author{
Alessandro Barbiero, Pier Alda Ferrari
}
\seealso{
\code{\link{ordcont}}, \code{\link{ordsample}}, \code{\link{corrcheck}}
}
\examples{
# Example 1
# consider a bivariate discrete random vector
k <- 2
# with these cumulative margins
marginal <- list(c(1/3, 2/3), c(0.1, 0.3, 0.6))
# generated discretizing a multivariate standard normal variable
# with correlation matrix
Sigma <- matrix(c(1, .75, .75, 1), 2, 2)
# the resulting joint distribution and correlation matrix
# for the bivariate discrete random vector are
res <- contord(marginal, Sigma, prob=TRUE)
res$pij
res$SigmaO
# let's check the margins are those assigned
cumsum(margin.table(res$pij,1))
cumsum(margin.table(res$pij,2))
# -> OK
# Example 2
# consider 4 discrete variables
k <- 4
# with these marginal distributions
marginal <- list(0.4,c(0.3,0.6), c(0.25,0.5,0.75), c(0.1,0.2,0.8,0.9))
# generated discretizing a multivariate standard normal variable
# with correlation matrix
Sigma <- matrix(0.5,4,4)
diag(Sigma) <- 1
# the resulting correlation matrix for the discrete variables is
contord(marginal, Sigma)
# note all the correlations are smaller than the original 0.6
# change Sigma, adding a negative correlation
Sigma[1,2] <- -0.15
Sigma[2,1] <- Sigma[1,2]
Sigma
# checking whether Sigma is still positive definite
eigen(Sigma)$values # all >0, OK
contord(marginal, Sigma)
# Example 2
# the same margins and the same correlation matrix
# but now we consider a 4-variate Students's t with df=3
contord(marginal, Sigma, df=3)
# -> a slight reduction in magnitude for all the correlations
}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
\keyword{distribution}
\keyword{datagen}
