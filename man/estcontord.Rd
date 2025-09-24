\name{estcontord}
\alias{estcontord}
\title{
Estimating based on a bivariate sample of the multivariate latent standard normal or Student's t distribution
}
\description{
Limited to the bivariate case, based on an iid sample, the function estimates the parameters of the t-copula-based discrete distribution, according to either a two-step approach (suggested) where the only unknown parameters are the copula correlation and the degrees of Freedom, whereas the marginal probabilities are estimated via the ecdf, or a full-maximum-likelihood approach, where also the two marginal probabilities are estimated jointly with the correlation and degrees of freedom
}
\usage{
estcontord(x, method="2-step")
}
\arguments{
  \item{x}{
a matrix with two columns containing integer values; it is the bivariate iid sample
}
  \item{method}{
\code{"2-step"} for estimating \code{rho} and \code{df} parametrically, and the margins via the ecdf; otherwise, a full-maximum-likelihood approach is carried out, where the marginal probabilities of the two random components are estimated as well
}
}
\value{
a list containing the estimates and for the "2-step" method their standard errors
}
\author{
Alessandro Barbiero, Pier Alda Ferrari
}
\seealso{
\code{\link{ordcont}},\code{\link{contord}}, \code{\link{ordsample}}, \code{\link{corrcheck}}
}
\examples{
## not run
#x1 <- c(rep(0,223),rep(1,269),rep(2,8))
#x2 <- c(rep(0,153), rep(1,70), rep(0,75),rep(1,187),
#        rep(2,7),rep(0,2),rep(1,4),rep(2,2))
#cor(x1,x2)
#x<-cbind(x1,x2)
#res <- estcontord(x)
#res
}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
\keyword{distribution}
\keyword{datagen}

