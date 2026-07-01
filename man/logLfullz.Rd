\name{logLfullz}
\alias{logLfullz}
\title{
Log-likelihood value for the t-copula-based model (with marginal probabilities transformed into unconstrained variables)
}
\description{
Negative log-likelihood value for the t-copula-based model to be used for estimation in the bivariate case
}
\usage{
logLfullz(arg, x)
}
\arguments{
  \item{arg}{
a vector containing the values of the correlation parameter \code{rho}, of the softmax-transformed \eqn{m_1} marginal probabilities for rv \eqn{X_1}, of the softmax-transformed \eqn{m_2} marginal probabilities for rv \eqn{X_2}, and, finally, of the degrees-of-freedom parameter \code{df}
}
  \item{x}{
a matrix with two columns containing the bivariate discrete data
}
}
\details{
The vector \code{arg} contains in position \eqn{1} the correlation parameter \code{rho}, from position \eqn{2} to position \eqn{m_1+1} the \eqn{m_1} values \eqn{z_i} from which the marginal probabilities of \eqn{X_1} are obtained as \eqn{p_{1i}=\exp(z_i)/\sum_{i=1}^{m_1} exp(z_i)}; from position \eqn{m_1+2} to position \eqn{m_1+m_2+1} the \eqn{m_2} values \eqn{w_j} from which the marginal probabilities of \eqn{X_2} are obtained as \eqn{p_{2j}=\exp(w_j)/\sum_{j=1}^{m_2} exp(w_j)}. The position \eqn{m_1+m_2+2} contains, if applicable, the degrees-of-freedom parameter.
}
\value{
the negative value of the log-likelihood function
}
\author{
Alessandro Barbiero, Pier Alda Ferrari
}
\seealso{
\code{\link{logLfull}}, \code{\link{estcontord}}
}
\examples{
x1 <- c(rep(0,223),rep(1,269),rep(2,8))
x2 <- c(rep(0,153), rep(1,70), rep(0,75),rep(1,187),
        rep(2,7),rep(0,2),rep(1,4),rep(2,2))
cor(x1,x2)
x <- cbind(x1,x2)
logLfullz(arg=c(0.5,rep(c(-0.5, 0, -3),2),5), x=x)
}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
\keyword{distribution}
\keyword{datagen}

