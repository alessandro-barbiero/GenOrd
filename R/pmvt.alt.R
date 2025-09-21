pmvt.alt <- function(low, upp, corr, df)
{
  integrand <- function(arg, df)
  {
    x <- arg[1]
    y <- arg[2]
    ff <- dmvt(c(x,y), sigma=corr, df=df, log=FALSE)
    return(ff)
  }
cuhre(f=integrand, df=df, lowerLimit=low, upperLimit=upp)$integral
}