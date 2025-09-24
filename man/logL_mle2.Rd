\name{logL_mle2}
\alias{logL_mle2}
\title{
Log-likelihood function for the t-copula-based model
}
\description{
Log-likelihood function for the t-copula-based model to be used for estimation in the bivariate case
}
\usage{
logL_mle2(rho, df, x)
}
\arguments{
  \item{rho}{
the value of \code{rho}
}
  \item{df}{
the value of \code{df}
}
  \item{x}{
a matrix with two columns containg the discrete data
}
}
\value{
the value of the log-likelihood function
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
x<-cbind(x1,x2)
logL_mle2(rho=0.5, df=5, x=x)
}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
\keyword{distribution}
\keyword{datagen}

