#' @importFrom stats qt
#' @importFrom utils head
#' @importFrom cubature cuhre
#' @importFrom mvtnorm pmvt
#' @importFrom mvtnorm dmvt
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats qnorm
#' @export
# (full) log-likelihood function for the t-copula-based discrete distribution
# function of rho, df and the (m1-1) plus (m2-1) probabilities of the two margins
logLfullz <- function(arg, x)
{
  k <- dim(x)[2]
  n <- dim(x)[1]
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
  p <- matrix(0,k1,k2)
  T <-table(x[,1],x[,2])
  rho <- arg[1]
  zx <- arg[2:(k1+1)]
  px <- exp(zx)/sum(exp(zx))
  zy <- arg[(k1+2):(k1+k2+1)]
  py <- exp(zy)/sum(exp(zy))
  Px <- c(0,cumsum(px));Px[k1+1]<-1
  Py <- c(0,cumsum(py));Py[k2+1]<-1
  Sigma <- matrix(c(1,rho,rho,1),2,2)
  if(length(arg)==k1+k2+2){
    df <- arg[k1+k2+2] # Student's t degrees of freedom
    if (df <= 0 | abs(rho)>1) return(1e06*n)
    thx <- qt(Px, df=df)
    thy <- qt(Py, df=df)
    for(i in 1:k1)
    {
      for(j in 1:k2)
      {
        # NEW FUNCTION
        P[i,j] <- p_rect_t(thx[i], thx[i+1], thy[j], thy[j + 1], rho, df)
      }
    }
    dl <- -sum(T*log(P))
  } else {
    thx <- qnorm(Px)
    thy <- qnorm(Py)
    for(i in 1:k1)
    {
      for(j in 1:k2)
      {
        P[i,j] <- pmvnorm(lower=c(thx[i],thy[j]), upper=c(thx[i+1],thy[j+1]), corr=Sigma)
      }
    }
    dl <- -sum(T*log(P))
  }
  return(dl)
}
