\name{estcontord}
\alias{estcontord}
\title{
Estimation based on a bivariate sample of the multivariate discrete model obtained by discretizing a latent standard normal or Student's t distribution
}
\description{
Limited to the bivariate case, based on an iid sample, the function estimates the parameters of the t-copula-based discrete distribution, according to either a two-step approach (suggested) where the only unknown parameters are the copula correlation and the degrees of freedom, whereas the two marginal probability distributions are estimated via the ecdf; or a full-maximum-likelihood approach, where also the two sets of marginal probabilities are estimated jointly with the copula correlation and degrees of freedom
}
\usage{
estcontord(x, method="2-step", start.df=5, control= list(), normal=FALSE)
}
\arguments{
  \item{x}{
a matrix with two columns containing integer values; it is the bivariate iid sample
}
  \item{method}{
\code{"2-step"} for estimating the Student's parameters \code{rho} and \code{df} parametrically, after the margins have been estimated via the univariate ecdf; otherwise, a full-maximum-likelihood approach is carried out, where the marginal probabilities of the two random components are estimated as well
}
 \item{start.df}{
starting value for the degrees-of-freedom parameter (by default, set equal to 5)
}
 \item{control}{
a list of control options for the \code{optim} function in the maximum likelihood optimization routine
}
\item{normal}{
a logical value; if \code{FALSE}, we estimate a discretized latent Student's t distribution; if \code{TRUE}, we estimate a discretized latent normal distribution
}
}
\value{
a list containing the estimates (and for the "2-step" method their standard errors) along with the maximum log-likelihood value
}
\author{
Alessandro Barbiero, Pier Alda Ferrari
}
\seealso{
\code{\link{logLfull}},\code{\link{logL_mle2}},\code{\link{ordcont}},\code{\link{contord}}, \code{\link{ordsample}}
}
\examples{
# EX.1
x1 <- c(rep(0,223),rep(1,269),rep(2,8))
x2 <- c(rep(0,153), rep(1,70), rep(0,75),rep(1,187),
        rep(2,7),rep(0,2),rep(1,4),rep(2,2))
cor(x1,x2)
x <- cbind(x1,x2)
res <- estcontord(x)
res
# EX.2
set.seed(12345)
rho <- 0.5
df <- 5
Sigma <- matrix(c(1,rho,rho,1), 2, 2)
margins <- list(c(1/3,2/3), c(1/6,1/2,5/6))
# drawing a sample of size n=500
# from a discretized bivariate t distribution with rho=0.5 and df=5
# with the assigned margins above
x <- ordsample(n=500, marginal=margins, Sigma=Sigma, cormat="continuous", df=df)
estcontord(x, "two-step")
# rho.hat ~ 0.477, df.hat ~ 7.02
estcontord(x, "full", control=list(maxit=2000, factr = 1e6))
# rho.hat ~ 0.477, df.hat ~ 6.87
}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
\keyword{distribution}
\keyword{datagen}

