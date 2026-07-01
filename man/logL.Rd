\name{logL}
\alias{logL}
\title{
Log-likelihood value for the t-copula-based model
}
\description{
Negative log-likelihood value for the t-copula–based model to be used for two-step estimation in the bivariate case. The marginal probabilities are set equal to the corresponding relative frequencies of the input data set.
}
\usage{
logL(arg, x)
}
\arguments{
  \item{arg}{
a vector containing the values of \code{rho} and \code{df}
}
  \item{x}{
a matrix with two columns containing the bivariate discrete data
}
}
\value{
the negative log-likelihood value
}
\author{
Alessandro Barbiero, Pier Alda Ferrari
}
\seealso{
\code{\link{estcontord}}
}
\examples{
x1 <- c(rep(0,223),rep(1,269),rep(2,8))
x2 <- c(rep(0,153), rep(1,70), rep(0,75),rep(1,187),
        rep(2,7),rep(0,2),rep(1,4),rep(2,2))
r <- cor(x1,x2)
x <- cbind(x1,x2)
logL(arg=c(r, 4), x=x)
# it is the same as
logL_mle2(rho=r, df=4, x=x)
}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
\keyword{distribution}
\keyword{datagen}

