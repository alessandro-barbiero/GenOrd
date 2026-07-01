#
# R code accompanying the paper
# Modeling correlated discrete data through Student’s t copula
# This file is available at: https://github.com/alessandro-barbiero/GenOrd
#

library(GenOrd)

# generate a symmetric U-shaped or reverse U-shaped distribution with k categories

gen_prob_seq <- function(k, type="uni") { # reverse U-shaped
  stopifnot(k > 2)
  if(type=="uni")
  {
    m <- floor(k / 2)
    left <- 2^seq(0, m)  # from 2^0 up to 2^m
    if (k %% 2 == 0) {
      # Even k: mirror all of left
      full_seq <- c(head(left,-1), rev(head(left,-1)))
    } else {
      # Odd k: mirror all except the middle term to avoid duplication
      full_seq <- c(left, rev(left[1:m]))
    }
  }
  else if(type=="bi") # U-shaped
  {
    m <- ceiling(k / 2)
    left <- 2^(seq(m-1, 0))
    full_seq <- c(left, rev(left))
    if (k %% 2 == 1) {
      full_seq <- full_seq[-m]
    }
  }
  prob_seq <- full_seq / sum(full_seq)
  return(prob_seq)
}
# examples
gen_prob_seq(5, "uni")
gen_prob_seq(5, "bi")

#
# Types of discrete distributions used in the numerical study of SECTION 3.3
# FIGURE 1
#
op<-par()
par(mfrow=c(1,4),mai=c(0.25,0.25,0.15,0.25),mgp=c(1,.5,0))
k <- 5
ylim=c(0,0.6)
F <- (1:(k-1))/k                             # UNIFORM
p <- c(F [1],diff(F),1-F[k-1])
barplot(p,ylim=ylim,names.arg=1:k)
F <- head(cumsum(gen_prob_seq(k,"uni")),-1)  # symmetric reversed-U shaped
p <- c(F [1],diff(F),1-F[k-1])
barplot(p,ylim=ylim,names.arg=1:k)
F <- cumsum((2^((k-1):1))/(2^k-1))           # geometric
p <- c(F [1],diff(F),1-F[k-1])
barplot(p,ylim=ylim,names.arg=1:k)
F <- head(cumsum(gen_prob_seq(k,"bi")),-1)   # symmetric bimodal U-shaped
p <- c(F [1],diff(F),1-F[k-1])
barplot(p,ylim=ylim,names.arg=1:k)
par(op)
#
# Total Variation distance between two probability tables
#
TVdist <- function(pij, qij)
{
  .5*sum(abs(pij-qij))
}

######################################################
################## final correlation #################
# for a bivariate rv with identical discrete margins #
########### obtained from discretization of ##########
######## a bivariate Student's t distribution ########
############ as a function of rho and df #############
########## for different types of margins ############
######################################################
# -> Tables 1 to 4, SECTION 3.3
#
marginal <- list()
rhoC <- (1:9)/10
cat <- c(3,4,5,6,7,8,9,10)
vdf <- c(3,10,20,Inf)
rhoD <- matrix(0,length(vdf)*length(rhoC),length(cat))

for(h in 1:length(rhoC))
{
  rho <- rhoC[h]
  for(i in 1:length(vdf)){
    df <- vdf[i]
    for(j in 1:length(cat))
    {
      k <- cat[j]
      # uniform margins
      marginal[[1]] <- (1:(k-1))/k
      marginal[[2]] <- (1:(k-1))/k
      # asym
      marginal[[1]] <- cumsum((2^((k-1):1))/(2^k-1))
      marginal[[2]] <- cumsum((2^((k-1):1))/(2^k-1))
      # sym uni
      marginal[[1]] <- head(cumsum(gen_prob_seq(k,"uni")),-1)
      marginal[[2]] <- head(cumsum(gen_prob_seq(k,"uni")),-1)
      # sym bim
      marginal[[1]] <- head(cumsum(gen_prob_seq(k,"bi")),-1)
      marginal[[2]] <- head(cumsum(gen_prob_seq(k,"bi")),-1)
      Sigma <- matrix(rho, 2, 2)
      diag(Sigma) <- 1
      res <- contord(marginal=marginal, Sigma=Sigma, df=df, prob=TRUE)
      SigmaD <- res$SigmaO
      rhoD[i+4*(h-1),j] <- SigmaD[2]
    }
  }
}
print(round(rhoD,4))
#

###############################################################
##################### Optimization problem ####################
# determine which is the marginal distribution that maximizes #
# the final correlation for assigned t-copula correlation and #
################ degrees of freedom parameters ################
###############################################################
# SECTION 3.4
library(nloptr)
postscript("maxttcorr.eps",width=12,height=6,horizontal=TRUE)
kvec   <- 2:7 # n.of categories
rhovec <- c(.1,.2,.3,.4,.5,.6,.7,.8,.9) # correlations
df <- 3 # degrees of freedom
par(mfrow=c(length(rhovec),length(kvec)), mai=c(0.15,0.4,0.15,0), oma=c(0,5,0,0))
set.seed(12345)
max.it <- 0
for(i in 1:length(rhovec))
{
  for(j in 1:length(kvec))
  {
    rho <- rhovec[i]
    k <- kvec[j]
    fn1 <- function(z){
      p <- exp(z)/sum(exp(z))
      F <- cumsum(p)[-length(p)] # common cdf
      print(F)
      Sigma <- matrix(c(1,rho,rho,1), 2, 2) # correlation matrix
      marginal <- list(F, F)
      s <- sign(rho)+(rho==0)
      -s*contord(marginal, Sigma, df=df)[2] # final correlation (with changed sign)
    }
    res <- nloptr(x0=rep(0,k), eval_f=fn1, opts=list(algorithm="NLOPT_LN_COBYLA",ftol_rel=1e-10,xtol_rel=1e-10, maxeval=1000))
    res$objective
    it.st <- 0
    while(res$status<1 |res$status>3)
    {
      res <- nloptr(x0=rep(0,k)+rnorm(k), eval_f=fn1, opts=list(algorithm="NLOPT_LN_COBYLA",ftol_rel=1e-4,xtol_rel=1e-4, maxeval=1000))
      it.st <- it.st + 1
    }
    max.it <- max(max.it,it.st)
    res$objective
    p <- exp(res$solution)/sum(exp(res$solution))
    p
    if(which.max(p)>k/2) p <- rev(p)
    barplot(p, ylim=c(0,0.9)) # plotting the "optimal" distribution
    text(k/2+1/2, 0.4, round(-res$objective,3),pos=3) # writing the maximum correlation
    if(j==1) mtext(substitute(list(rho^{(t)}) == list(x),list(x = rho)),side=2,line=3,las=1)
  }
}
dev.off()


#
# SECTION 5.1
rho <- 0.5
df <- 3
Sigma <-matrix(c(1, rho, rho, 1), 2, 2)
F <- c(0.2, 0.4, 0.6, 0.8)
marginal <- list(F, F)
P.G <- contord(marginal, Sigma, prob=TRUE)
lapply(P.G, function(m) round(m, 4))
P.t <- contord(marginal, Sigma, df=3, prob=TRUE)
lapply(P.t, function(m) round(m, 4))
round(prop.table(P.G[[1]], margin = 1), 4)
round(prop.table(P.t[[1]], margin = 1), 4)

# SECTION 5.2
set.seed(12345)
k <- 7
F <- cumsum((2^((k-1):1))/(2^k-1))
marginal <- list(F,F)
corrcheck(marginal)
rho <- 0.5
Sigma <-matrix(c(1,rho,rho,1), 2, 2)
df <- 5
res.ord <- contord(marginal, Sigma, df=df, prob=TRUE)
round(res.ord$pij, 4)
res.ord$SigmaOrd
x <- ordsample(10000, marginal, Sigma, df=df, cormat="continuous")
head(x)
cor(x)
rho.ord <- res.ord$SigmaOrd[1,2]
ordcont(marginal, Sigma=matrix(c(1,rho.ord,rho.ord,1), 2, 2), df=df)
#


# SECTION 5.3
# simulated data set from
# K. F. Sellers, D. S. Morris, and N. Balakrishnan
# Bivariate Conway–Maxwell–Poisson distribution: formulation, properties, and inference
# Journal of Multivariate Analysis, 150:152–168, 2016.

x1 <- c(rep(0,223),rep(1,269),rep(2,8))
x2 <- c(rep(0,153), rep(1,70), rep(0,75),rep(1,187),
        rep(2,7),rep(0,2),rep(1,4),rep(2,2))
cor(x1,x2)
x<-cbind(x1,x2)
res <- estcontord(x) # two-step approach
res
rho.est <- res$estimates[1]
Sigma <- matrix(c(1,rho.est,rho.est,1),2,2)
marginal <- list(cumsum(res$margin1[-3]),cumsum(res$margin2[-3]))
res.check <- contord(marginal=marginal, Sigma, df=res$estimates[2], integerdf=FALSE, prob=TRUE)
res.check$pij # fitted joint distribution
prop.table(table(x[,1],x[,2])) # empirical joint distribution
TVdist(res.check$pij, prop.table(table(x[,1],x[,2])))

# SECTION 5.4
##########################
## Analysis of real data #
### Type D personality ###
##########################
# data set
# taken from
# L. Kolbe, F. Oort, and S. Jak
# Bivariate distributions underlying responses to ordinal variables
# Psych, 3(4):562–578, 2021

freq <- matrix(c(67,41,34,35,24,
                 15,28,48,22,10,
                 16,30,39,34,11,
                 8,4,11,28,11,
                 3,2,1,5,9),5,5)
n <- sum(freq)
tab <- as.table(freq)
X1 <- rep(1:5, times = rowSums(tab))
X2 <- unlist(lapply(1:5, function(i) rep(1:5, times = tab[i, ])))
X <- cbind(X1, X2)
cor(X)
# two-step approach
res.mdpi <- estcontord(X, method="2-step")
res.mdpi
margin1 <- cumsum(res.mdpi$margin1)[-5]
margin2 <- cumsum(res.mdpi$margin2)[-5]
rho <- res.mdpi$estimates[1]
df  <- res.mdpi$estimates[2]
# fitted joint distribution:
tab.est <- contord(list(margin1, margin2), matrix(c(1,rho,rho,1),2,2), df=df, integerdf=FALSE, prob=TRUE)$pij*n
# chi-square statistics
Chi <- sum((tab-tab.est)^2/tab.est)
Chi
1-pchisq(Chi, 5*5-1-2-2*4)
# likelihood ratio chi-square statistic
G2 <- 2*sum(tab*log(tab/tab.est))
G2
# full maximum likelihood approach
res.mdpi.full <- estcontord(X, method="full")
marginal <- list(cumsum(res.mdpi.full[3:7][-5]), cumsum(res.mdpi.full[8:12][-5]))
rho <- res.mdpi.full[1]
Sigma <- matrix(c(1,rho,rho,1),2,2)
df <- res.mdpi.full[2]
tab.est.full <- contord(marginal, Sigma, df=df, integerdf=FALSE, prob=TRUE)$pij*n
arg.full <- res.mdpi.full[c(1,3:12,2)]
arg.full[2:11]<-log(arg.full[2:11]) # invertirg the softmax function
logLfull(arg.full,x=X)
Chi.full <- sum((tab-tab.est.full)^2/tab.est.full)
Chi.full
G2.full <- 2*sum(tab*log(tab/tab.est.full))
G2.full
1-pchisq(Chi.full, 5*5-1-2-2*4)

