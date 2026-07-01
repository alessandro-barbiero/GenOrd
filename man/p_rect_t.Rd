\name{p_rect_t}
\alias{p_rect_t}
\title{Rectangle Probability for a Bivariate Student's t Distribution}
\description{
Computes the probability that a bivariate Student's t random vector with
correlation parameter \code{rho} and \code{df} degrees of freedom (non necessarily integer) falls in a
rectangular region.
}
\usage{
p_rect_t(ax, bx, ay, by, rho, df, rel.tol = 1e-4)
}
\arguments{
  \item{ax}{lower bound for the first component}
  \item{bx}{upper bound for the first component}
  \item{ay}{lower bound for the second component}
  \item{by}{upper bound for the second component}
  \item{rho}{correlation parameter of the bivariate Student's t distribution.}
  \item{df}{degrees-of-freedom parameter}
  \item{rel.tol}{relative accuracy requested for numerical integration; see
  \code{\link[stats]{integrate}}.}
}
\details{
The function evaluates \deqn{P(a_x < X \le b_x, a_y < Y \le b_y),} where \eqn{(X,Y)} follows a standard bivariate Student's t distribution with correlation parameter \code{rho} and \code{df} degrees of freedom.

The probability is computed through one-dimensional numerical integration after  transformation to the marginal t probability scale. Infinite lower and upper ounds are allowed.
}
\value{
A numeric value giving the rectangle probability.
}
\author{
Alessandro Barbiero
}
\seealso{
\code{\link{pmvt.alt}}, \code{\link[mvtnorm]{pmvt}},
\code{\link[stats]{integrate}}
}
\examples{
## Probability over a finite rectangle
p_rect_t(ax = -1, bx = 0, ay = -0.5, by = 1, rho = 0.5, df = 3.5)

## Probability with an infinite lower bound
p_rect_t(ax = -Inf, bx = 0, ay = -0.5, by = 1, rho = 0.5, df = 5)
}
\keyword{distribution}
\keyword{multivariate}
