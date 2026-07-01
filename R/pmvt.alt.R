#' @export
pmvt.alt <- function (low, upp, corr, df)
{
  integrand <- function(arg, df) {
    ff <- dmvt(arg, sigma = corr, df = df, log = FALSE)
    return(ff)
  }
  cuhre(f = integrand, df = df, lowerLimit = low, upperLimit = upp)$integral
}
