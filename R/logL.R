#' @importFrom stats qt
#' @importFrom utils head
#' @importFrom cubature cuhre
#' @importFrom mvtnorm dmvt
#' @export
#########################################################################
### ESTIMATING a GAUSSIAN (or t) COPULA for (BIVARIATE) DISCRETE DATA ###
#########################################################################

# loglikelihood function the t-copula-based discrete distribution
# considering only the copula parameters rho and df as the unknown parameters
logL <- function(arg, x)
{
  rho <- arg[1]
  df  <- arg[2]

  k <- dim(x)[2]
  n <- dim(x)[1]
  if (df <= 0 || abs(rho) >= 1) return(n*1e6)

  emarginal <-  vector("list", k)
  esupport <-  vector("list", k)

  for(i in 1:k)
  {
    n <- length(x[,1])
    t<-table(x[,i])
    emarginal[[i]]<-head(cumsum(as.numeric(t))/n,-1)
    esupport[[i]]<-as.numeric(attr(t,"dimnames")[[1]])
  }

  k1<-length(esupport[[1]])
  k2<-length(esupport[[2]])
  P <- matrix(0,k1,k2)
  T <-table(x[,1],x[,2])

  thx <- c(-Inf,qt(emarginal[[1]],df),Inf)
  thy <- c(-Inf,qt(emarginal[[2]],df),Inf)
  if (df <= 0) return(999*n)
  Sigma <- matrix(c(1,rho,rho,1),2,2)

  integrand <- function(arg, df)
  {
    x <- arg[1]
    y <- arg[2]
    ff <- dmvt(c(x,y), sigma=Sigma, df=df, log=FALSE)
    return(ff)
  }

  for(i in 1:k1)
  {
    for(j in 1:k2)
    {
      P[i,j] = cuhre(f=integrand, df=df, lowerLimit=c(thx[i],thy[j]), upperLimit=c(thx[i+1],thy[j+1]))$integral
      #pStudent(lower=c(thx[i],thy[j]),upper=c(thx[i+1],thy[j+1]),df=df,loc=c(0,0),scale=Sigma)
    }
  }
  if (any(P <= 0)) return(n*1e6)

  dl <- -sum(T*log(P))
  return(dl)
}
