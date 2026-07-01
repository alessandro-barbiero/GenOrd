\name{pmvt.alt}
\alias{pmvt.alt}
\title{
Probabilities for a multivariate t with non-integer degrees-of-freedom parameter
}
\description{
Computes probabilities of an hyper-rectangle for a standard multivariate Student's \eqn{t} with non-integer degrees-of-freedom parameter, directly integrating the function \code{dmvt}, available in the package \code{mvtnorm}, through the \code{cuhre} function, available in the package \code{cubature}
}
\usage{
pmvt.alt(low, upp, corr, df)
}
\arguments{
  \item{low}{
a vector containing the lower bounds
}
  \item{upp}{
a vector containing the upper bounds
}
\item{corr}{
the correlation matrix of the \eqn{t} copula
}
\item{df}{
the degrees-of-freedom parameter, possibly non-integer
}
}
\author{
Alessandro Barbiero, Pier Alda Ferrari
}
\seealso{
\code{\link{contord}}
}
\examples{
# dimension 2
lower <- c(-1, -1/2)
upper <- c(0, 1)
Sigma <- matrix(c(1, 0.5, 0.5, 1), 2 , 2)
pmvt.alt(low=lower, upp=upper, corr=Sigma, df=3.5)
# dimension 4
lower <- c(-1, -1/2, 1/2, -1/4)
upper <- c(0, 1, 2, 3/4)
Sigma <- matrix(0.25, 4, 4)
diag(Sigma) <- 1
pmvt.alt(low=lower, upp=upper, corr=Sigma, df=3.5)
# let's set df=Inf
pmvt.alt(low=lower, upp=upper, corr=Sigma, df=Inf)
# this probability should be equal to
mvtnorm::pmvnorm(lower=lower, upper=upper, corr=Sigma) # OK!
}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
\keyword{distribution}
\keyword{datagen}

