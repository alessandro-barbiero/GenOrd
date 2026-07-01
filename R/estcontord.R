#' @importFrom stats cor optim
#' @importFrom cubature cuhre
#' @importFrom bbmle mle2
#' @importFrom mvtnorm dmvt
#' @export
#########################################################################
### ESTIMATING a GAUSSIAN (or t) COPULA for (BIVARIATE) DISCRETE DATA ###
#########################################################################

# Two-Step or Full Maximum Likelihood Estimation for t-copula-based bivariate discrete data
#
estcontord <- function(x, method="2-step", start.df=5, control=list(),  normal=FALSE)
{
  n <- dim(x)[1]
  marg1 <- table(x[,1])/n
  marg2 <- table(x[,2])/n
  if(method=="2-step"||method=="two-step")
  {
    res <- mle2(minuslogl = logL_mle2, start = list(rho = cor(x)[1,2], df = start.df), data = list(x = x), control=control)
    ret <- list(margin1=marg1, margin2=marg2, estimates=c(res@coef,se=sqrt(diag(res@vcov)),log.lik.max=res@min))
    return(ret)
  }
  else
  {
    emarginal <-  vector("list", 2)
    esupport  <-  vector("list", 2)
    ez <- vector("list",2)
    for(i in 1:2)
    {
      t <- table(x[,i])
      emarginal[[i]] <- prop.table(table(x[,i]))#head(cumsum(as.numeric(t))/n,-1)
      esupport[[i]]  <- as.numeric(attr(t,"dimnames")[[1]])
      ez[[i]] <- log(emarginal[[i]]/emarginal[[i]][1])
    }
    k1 <- length(esupport[[1]])
    k2 <- length(esupport[[2]])
    #
    if(normal==FALSE){
      par <- c(cor(x)[2],ez[[1]],ez[[2]], start.df)
      #lower <- c(-.999,rep(-Inf,length(par)-2),1)
      #upper <- c(.999,rep(Inf,length(par)-2),Inf)
    } else {
      par <- c(cor(x)[2],ez[[1]],ez[[2]])
      #lower <- c(-.999,rep(-Inf,length(par)-1))
      #upper <- c(.999,rep(Inf,length(par)-1))
    }
    #res <- optim(par=par, fn=logLfullz, x=x, lower=lower, upper=upper, method="L-BFGS-B", control = list(factr = 1e6))
    res <- optim(par=par, fn=logLfullz, x=x, control=control)
    z1 <- res$par[2:(k1+1)]
    z2 <- res$par[(k1+2):(k1+k2+1)]
    p1 <- exp(z1)/sum(exp(z1))
    p2 <- exp(z2)/sum(exp(z2))
    return(c(rho=res$par[1], df=res$par[k1+k2+2], margin1=p1, margin2=p2, log.lik.max=res$value))
  }
}
