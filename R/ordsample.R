#' @export
#' @importFrom stats qt
#' @importFrom mvtnorm rmvt
#' @importFrom mvtnorm rmvnorm
ordsample <-
function(n, marginal, Sigma, support=list(), Spearman=FALSE, cormat="discrete", df=Inf)
{
    if (!all(unlist(lapply(marginal, function(x) (sort(x)==x & min(x)>0 & max(x)<1))))) stop("Error in assigning marginal distributions!")
    if(!isSymmetric(Sigma) | min(eigen(Sigma)$values)<0 | !all(diag(Sigma)==1)) stop("Correlation matrix not valid!")
# k=number of variables
k <- length(marginal)
# kj=number of categories for the k variables (vector of k integer numbers)
kj <- numeric(k)
len <- length(support)
for(i in 1:k)
{
kj[i] <- length(marginal[[i]]) + 1
if(len==0)
{
support[[i]] <- 1:kj[i]
}
}
if(cormat=="discrete")
{
Sigmac <- ordcont(marginal=marginal, Sigma=Sigma, support=support, Spearman=Spearman, df=df)[[1]]
Sigma <- Sigmac
}
# sample of size n from k-dimensional Student's t or normal with vector of zero means and correlation matrix Sigma
valori <- rmvt(n, sigma=Sigma, delta=rep(0, k), df=df, type="shifted")
if(n==1) valori <- matrix(valori,nrow=1)
# discretization according to the marginal distributions
for(i in 1:k)
{
valori[,i] <- as.integer(cut(valori[,i], breaks=c(min(valori[,i])-1, qt(marginal[[i]], df=df), max(valori[,i])+1)))
valori[,i] <- support[[i]][valori[,i]]
}
return(valori)
}

