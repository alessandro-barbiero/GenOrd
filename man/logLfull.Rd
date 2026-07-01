\name{logLfull}
\alias{logLfull}
\title{
Log-likelihood value for the t-copula-based model
}
\description{
Negative log-likelihood value for the t-copula-based model, to be used for estimation in the bivariate case
}
\usage{
logLfull(arg, x)
}
\arguments{
  \item{arg}{
a vector containing the values of the correlation parameter \code{rho}, of the \eqn{m_1} marginal probabilities for rv \eqn{X_1}, of the \eqn{m_2} marginal probabilities for rv \eqn{X_2}, and, finally, of the degrees-of-freedom parameter \code{df}
}
  \item{x}{
a matrix with two columns containing the bivariate discrete data
}
}
\details{
The vector \code{arg} contains in position \eqn{1} the correlation parameter \code{rho}, from position \eqn{2} to position \eqn{m_1+1} the \eqn{m_1} values of the marginal probabilities of \eqn{X_1} ; from position \eqn{m_1+2} to position \eqn{m_1+m_2+1} the \eqn{m_2} values of the marginal probabilities of \eqn{X_2}. The position \eqn{m_1+m_2+2} contains, if applicable, the degrees-of-freedom parameter.
}
\value{
the negative value of the log-likelihood function
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
cor(x1,x2)
x <- cbind(x1,x2)
logLfull(arg=c(0.5,rep(c(0.44, 0.54, 0.02),2),5), x=x)
}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
\keyword{distribution}
\keyword{datagen}

