\name{pmvt.alt}
\alias{pmvt.alt}
\title{
Probabilities for a multivariate t with non-integer degrees-of-freedom parameter
}
\description{
Computes probabilities for a multivariate \eqn{t} with non-integer degrees-of-freedom parameter, directly integrating the function \code{dmvt}
}
\usage{
pmvt.alt(low, upp, corr, df)
}
\arguments{
  \item{low}{
a vector containing the bivariate lower bounds
}
  \item{upp}{
a vector containing the bivariate upper bounds
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
margin1 <- c(0.2,0.4,0.6,0.8)
margin2 <- margin1
marginal <- list(margin1, margin2)
sigma <- matrix(c(1,0.3,0.3,1),2,2)
df <- 3.5
contord(marginal=marginal, Sigma=sigma, df=df, integerdf=FALSE, prob=TRUE)
# compare with
contord(marginal=marginal, Sigma=sigma, df=round(df), prob=TRUE)
}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
\keyword{distribution}
\keyword{datagen}

