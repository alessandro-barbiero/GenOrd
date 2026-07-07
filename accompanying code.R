################################################################################
# R code accompanying the paper                                                #
# Modeling correlated discrete data                                            #
# through the multivariate Student's t distribution                            #
################################################################################

library(GenOrd)

###########################
# some numerical examples #
###########################
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

#
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

# degrees-of-freedom and tail dependence
# in the discrete setting
rho <- 0.5
Sigma <- matrix(c(1,rho,rho,1),2,2)
m <- 4
# discrete uniform margins
margins <- list((1:3)/4, (1:3)/4)
contord(margins, Sigma, df=1, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=3, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=10, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=20, prob=TRUE)$pij[m,m]
# reverse-U shaped margins
margins <- list(c(1,3,5)/6,c(1,3,5)/6)
contord(margins, Sigma, df=1, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=3, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=10, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=20, prob=TRUE)$pij[m,m]
# U shaped margins
margins <- list(c(2,3,4)/6,c(2,3,4)/6)
contord(margins, Sigma, df=1, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=3, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=10, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=20, prob=TRUE)$pij[m,m]
# margins with decreasing pattern
F <- cumsum((2^((4-1):1))/(2^4-1))
margins <- list(F, F)
contord(margins, Sigma, df=1, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=3, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=10, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=20, prob=TRUE)$pij[m,m]
# Increasing the degrees-of-freedom parameter,
# while keeping the margins and the latent correlation fixed,
# decreases the probability of the (m,m) cell,
# which may be viewed as a rough measure of tail dependence
# in the discrete setting

# also with unequal margins
margins <- list((1:3)/4, c(1,3,5)/6)
contord(margins, Sigma, df=1, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=3, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=10, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=20, prob=TRUE)$pij[m,m]

margins <- list((1:3)/4, c(2,3,4)/6)
contord(margins, Sigma, df=1, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=3, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=10, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=20, prob=TRUE)$pij[m,m]

margins <- list((1:3)/4, cumsum((2^((4-1):1))/(2^4-1)))
contord(margins, Sigma, df=1, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=3, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=10, prob=TRUE)$pij[m,m]
contord(margins, Sigma, df=20, prob=TRUE)$pij[m,m]
# etc.

###########################
# Example for Section 2.2 #
###########################
# uniform margins with 5 categories
margins <- list(c(1,2,3,4)/5,c(1,2,3,4)/5)
Sigma <- matrix(c(1,0.5,0.5,1),2,2)
df.1 <- 20
df.2 <- 2
# 
Sigma.1 <- ordcont(margins, Sigma, df=df.1)$SigmaC
rho.1 <- Sigma.1[1,2]
rho.1
Sigma.2 <- ordcont(margins, Sigma, df=df.2)$SigmaC
rho.2 <- Sigma.2[1,2]
rho.2
# check
contord(margins, Sigma.1, df=df.1)
contord(margins, Sigma.2, df=df.2)
# rho^t.1 and df.1 yields the same correlation as rho^t.1 and df.2


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
# Types of discrete distributions used in the numerical study of SECTION 2.3
# FIGURE 1
#
op<-par()
par(mfrow=c(1,4),mai=c(0.25,0.25,0.15,0.25),mgp=c(1,.5,0))
k <- 5
ylim=c(0,0.6)
F <- (1:(k-1))/k                             # UNIFORM
p <- c(F [1],diff(F),1-F[k-1])
barplot(p,ylim=ylim,names.arg=1:k)
F <- cumsum((2^((k-1):1))/(2^k-1))           # DECREASING (geometric)
p <- c(F [1],diff(F),1-F[k-1])
barplot(p,ylim=ylim,names.arg=1:k)
F <- head(cumsum(gen_prob_seq(k,"uni")),-1)  # symmetric REVERSED-U-SHAPED
p <- c(F [1],diff(F),1-F[k-1])
barplot(p,ylim=ylim,names.arg=1:k)
F <- head(cumsum(gen_prob_seq(k,"bi")),-1)   # symmetric (anti-unimodal) U-SHAPED
p <- c(F [1],diff(F),1-F[k-1])
barplot(p,ylim=ylim,names.arg=1:k)
par(op)
#
# Distances between Probability Tables
# Total Variation distance
TVdist <- function(pij, qij)
{
  .5*sum(abs(pij-qij))
}
# Hellinger distance
Helldist <- function(pij, qij)
{
  sqrt(.5*sum((sqrt(pij)-sqrt(qij))^2))
}
# Euclidean distance
Eudist <- function(pij, qij)
{
  sqrt(sum((pij-qij)^2))
}

######################################################
################## final correlation #################
# for a bivariate rv with identical discrete margins #
########### obtained from discretization of ##########
######## a bivariate Student's t distribution ########
############ as a function of rho and df #############
########## for different types of margins ############
######################################################
# -> Section 2.3, Table 1, and
# -> Tables B1 to B4, appendix B
#
marginal <- list()
rhoC <- (1:9)/10
cat <- c(3,4,5,6,7,8,9,10) # add 2 for uniform and asymmetrical margins!
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
# -> SECTION 2.3, Table 2, and Appendix C, Figure C1 to C4
library(nloptr)
# postscript("maxttcorr.eps",width=12,height=6,horizontal=TRUE)
kvec   <- 2:7 # n.of categories
rhovec <- c(.1,.2,.3,.4,.5,.6,.7,.8,.9) # correlations

df <- Inf # degrees of freedom

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
      res <- nloptr(x0=rep(0,k)+rnorm(k), eval_f=fn1, opts=list(algorithm="NLOPT_LN_COBYLA",
                                    ftol_rel=1e-4,xtol_rel=1e-4, maxeval=1000))
      it.st <- it.st + 1
    }
    max.it <- max(max.it,it.st)
    res$objective
    p <- exp(res$solution)/sum(exp(res$solution))
    p
    if(which.max(p)>k/2) p <- rev(p)
    barplot(p, ylim=c(0,0.9)) # plotting the "optimal" distribution
    text(k/2+1/2, 0.4, formatC(-res$objective, format = "f", digits = 3), pos=3, cex=1.25)
    # writing the maximum correlation
    if(j==1) mtext(substitute(list(rho^{(t)}) == list(x),list(x = rho)),side=2,line=3,las=1)
  }
}
# dev.off()


# -> SECTION 2.4, Table 3
# evaluating Total Variation Distance between Gaussian-copula and t-copula-based
# joint distributions, as functions of nu, rho, m and type of margin

# uniform margins
marginal <- list()
rhoC <- (1:9)/10
cat <- c(2,3,4,5,6,7,8,9,10)
vdf <- c(3,10,20,Inf)
TV <- matrix(0,length(vdf)*length(rhoC),length(cat))
for(h in 1:length(rhoC))
{
  rho <- rhoC[h]
  for(i in 1:length(vdf)){
    df <- vdf[i]
    for(j in 1:length(cat))
    {
      k <- cat[j]
      # uniform margins (or any other margin...)
      marginal[[1]] <- (1:(k-1))/k
      marginal[[2]] <- (1:(k-1))/k
      Sigma <- matrix(rho, 2, 2)
      diag(Sigma) <- 1
      res <- contord(marginal=marginal, Sigma=Sigma, df=df, prob=TRUE)
      resinf <- contord(marginal=marginal, Sigma=Sigma, df=Inf, prob=TRUE)
      SigmaD <- res$SigmaO
      TV[i+4*(h-1),j] <- TVdist(res$pij,resinf$pij)
    }
  }
}
print(round(TV[-4*(1:9),],4))


##########################
# SECTION 4              #
# Data analyses          #
#                        #
##########################
####### SECTION 4.1 ######
##########################
# simulated data set from
# K. F. Sellers, D. S. Morris, and N. Balakrishnan
# Bivariate Conway–Maxwell–Poisson distribution: formulation, properties, and inference
# Journal of Multivariate Analysis, 150:152–168, 2016.

x1 <- c(rep(0,223),rep(1,269),rep(2,8))
x2 <- c(rep(0,153), rep(1,70), rep(0,75),rep(1,187),
        rep(2,7),rep(0,2),rep(1,4),rep(2,2))
cor(x1,x2)
x <- cbind(x1,x2)
tab <- table(x1,x2)
res <- estcontord(x) # two-step approach
res

rho.est <- res$estimates[1]
Sigma <- matrix(c(1,rho.est,rho.est,1),2,2)
marginal <- list(cumsum(res$margin1[-3]),cumsum(res$margin2[-3]))
res.check <- contord(marginal=marginal, Sigma, df=res$estimates[2], integerdf=FALSE, prob=TRUE)
res.check$pij # fitted joint distribution
tab.obs <- prop.table(table(x[,1],x[,2])) # empirical joint distribution
# distances between probability tables
TVdist(res.check$pij, tab.obs)
Helldist(res.check$pij, tab.obs)
Eudist(res.check$pij, tab.obs)

res.full <- estcontord(x, method="full") # full-likelihood maximization
res.full
margin1 <- cumsum(res.full[3:5])[-3]
margin2 <- cumsum(res.full[6:8])[-3]
rho <- res.full[1]
df  <- res.full[2]
# fitted joint distribution:
tab.est.full <- contord(list(margin1, margin2), matrix(c(1,rho,rho,1),2,2), df=df,
                        integerdf=FALSE, prob=TRUE)$pij
round(tab.est.full, 4)
TVdist(tab.est.full, tab.obs)
# fitting the discretized Gaussian model
res.G <- estcontord(x, method="full", normal=TRUE)
res.G
rho <- res.G[1]
Sigma <- matrix(c(1,rho,rho,1),2,2)
margin1 <- cumsum(res.G[3:4])
margin2 <- cumsum(res.G[6:7])
tab.G   <- contord(list(margin1,margin2), Sigma, prob=TRUE)$pij
TVdist(tab.G, tab.obs)

##########################
####### SECTION 4.2 ######
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
# fitting the discretized Gaussian
res.G <- estcontord(X, method="full", normal=TRUE, control=list(maxit=5000))
res.G
rho <- res.G[1]
Sigma <- matrix(c(1,rho,rho,1),2,2)
margin1 <- cumsum(res.G[3:6])
margin2 <- cumsum(res.G[8:11])
tab.G   <- contord(list(margin1,margin2), Sigma, prob=TRUE)$pij*n
tab.G
sum((tab-tab.G)^2/tab.G)
TVdist(tab.G/n, tab/n)
# two-step approach
res.mdpi <- estcontord(X, method="2-step")
res.mdpi
margin1 <- cumsum(res.mdpi$margin1)[-5]
margin2 <- cumsum(res.mdpi$margin2)[-5]
rho <- res.mdpi$estimates[1]
df  <- res.mdpi$estimates[2]
# fitted joint distribution:
tab.est <- contord(list(margin1, margin2), matrix(c(1,rho,rho,1),2,2), df=df,
                   integerdf=FALSE, prob=TRUE)$pij*n
# chi-square statistic
Chi <- sum((tab-tab.est)^2/tab.est)
Chi
1-pchisq(Chi, 5*5-1-2-2*4)
# likelihood ratio chi-square statistic
G2 <- 2*sum(tab*log(tab/tab.est))
G2
# full maximum likelihood approach
res.mdpi.full <- estcontord(X, method="full", control=list(maxit=5000))
res.mdpi.full
marginal <- list(cumsum(res.mdpi.full[3:7][-5]), cumsum(res.mdpi.full[8:12][-5]))
rho <- res.mdpi.full[1]
Sigma <- matrix(c(1,rho,rho,1),2,2)
df <- res.mdpi.full[2]
tab.est.full <- contord(marginal, Sigma, df=df, integerdf=FALSE, prob=TRUE)$pij*n
Chi.full <- sum((tab-tab.est.full)^2/tab.est.full)
Chi.full
G2.full <- 2*sum(tab*log(tab/tab.est.full))
G2.full
1-pchisq(Chi.full, 5*5-1-2-2*4)
# TV distance
TVdist(tab.est/n, tab/n)

##########################
####### SECTION 4.3 ######
##########################
#  Analysis of real data #
#     Premarital sex     #
#          and           #
#  Teenage birth control #
##########################
# Example from Agresti,
# An Introduction to Categorical Data Analysis (2018)
# table 7.10, p.214
tab <- matrix(c(81,24,18,36,68,26,41,57,60,29,74,161,38,14,42,157),4,4)
n <- sum(tab)
x <- matrix(NA, n, 2)
# building the bivariate sample as an n times 2 matrix
k <- 1
for (i in 1:nrow(tab)) {
  for (j in 1:ncol(tab)) {
    f <- tab[i, j]
    if (f > 0) {
      x[k:(k+f-1), ] <- cbind(rep(i, f), rep(j, f))
      k <- k + f
    }
  }
}
res.agresti <- estcontord(x, "full", control=list(maxit=5000))
res.agresti
rho     <- res.agresti[1]
Sigma   <-matrix(c(1,rho,rho,1),2,2)
df      <-  res.agresti[2]
margin1 <- cumsum(res.agresti[3:5])
margin2 <- cumsum(res.agresti[7:9])
tab.est <- GenOrd:::contord(list(margin1,margin2), Sigma, df=df, integerdf=FALSE, prob=TRUE)$pij*n
tab.est
chi.t <- sum((tab.est-tab)^2/tab.est)
chi.t
1-pchisq(chi.t, 4*4-4-4-1)
TVdist(tab/n, tab.est/n)
# discretized Gaussian model
res.G <- estcontord(x, "full", normal=TRUE, control=list(maxit=5000))
res.G
rho <- res.G[1]
Sigma <- matrix(c(1,rho,rho,1),2,2)
margin1 <- cumsum(res.G[3:5])
margin2 <- cumsum(res.G[7:9])
tab.n   <- contord(list(margin1,margin2), Sigma, prob=TRUE)$pij*dim(x)[1]
tab.n
chi.n   <- sum((tab-tab.n)^2/tab.n)
1-pchisq(chi.n, 4*4-4-4)
TVdist(tab/n, tab.n/n)

###############################
##         SECTION 3         ##
# additional simulation study #
## for assessing statistical ##
## properties of estimators  ##
###############################
margin.1 <- c(0.25, 0.5, 0.75)            # Uniform
margin.2 <- cumsum((2^((4-1):1))/(2^4-1)) # Decreasing
margin.2 <- c(1,3,5)/6                    # reverse-U shaped
margin.2 <- c(2,3,4)/6                    # U-shaped
# Uniform - Decreasing - Reverse U-shaped - U-shaped
margin.type <- "U-U"
margin.2 <- margin.1
# t correlation (-0.75, -0.25, 0.25, 0.75)
rho <- 0.75
Sigma <- matrix(c(1,rho,rho,1),2,2)
# degrees of freedom (2, 5, 10)
nu <- 10
N <- 1000
n <- 100
filename <- paste("Latent_T_margins",margin.type,"n",n,"rho",rho,"nu",nu,"N",N)
est.nu <- numeric(N)
est.rho<- numeric(N)
set.seed(12345)
t1 <- Sys.time()
# TWO-STEP ESTIMATION
for(i in 1:N)
{
  x <- ordsample(n=n, marginal=list(margin.1,margin.2), Sigma=Sigma, df=nu, cormat="continuous")
  res<-estcontord(x, method="2-step", control=list(maxit=2000))
  est.rho[i]<-res$estimates[1]
  est.nu[i] <-res$estimates[2]
  print(i)
}
t2 <- Sys.time()
#
est.nu.full  <- numeric(N)
est.rho.full <- numeric(N)
set.seed(12345)
t3 <- Sys.time()
# FULL ML ESTIMATION
for(i in 1:N)
{
  x <- ordsample(n=n, marginal=list(margin.1,margin.2), Sigma=Sigma, df=nu, cormat="continuous")
  res<-tryCatch(estcontord(x, method="full", control=list(maxit=2000)),error = function(e) NaN)
  if(any(!is.finite(res))){
    # res<-tryCatch(estcontord(x, method="full", control=list(maxit=1000)),error = function(e) NaN)
    # if(any(!is.finite(res)))
    # {
      est.rho.full[i]<-NA
      est.nu.full[i] <-NA
     next
    #}
  }
  est.rho.full[i]<-res["rho"]
  est.nu.full[i] <-res["df"]
  print(i)
}
t4 <- Sys.time()
sink(file=paste(filename,"txt",sep="."))
cat("Margin type:\n")
margin.type
cat("Margins:\n")
margin.1
margin.2
cat("Mean rho.hat:\n")
mean(est.rho)
cat("Median rho.hat:\n")
median(est.rho)
cat("sd rho.hat:\n")
sd(est.rho)
cat("Quartiles nu.hat:\n")
quantile(est.nu, c(0.25,0.5,0.75))
#
cat("Mean rho.hat FULL:\n")
mean(est.rho.full, na.rm=TRUE)
cat("Median rho.hat FULL:\n")
median(est.rho.full, na.rm=TRUE)
cat("sd rho.hat FULL:\n")
sd(est.rho.full, na.rm=TRUE)
cat("Quartiles nu.hat FULL:\n")
quantile(est.nu.full, c(0.25,0.5,0.75), na.rm=TRUE)
# Time 2-step:
cat("Time (secs) per simulation 2-step: \n")
as.numeric(t2-t1,"secs")/N
# Time FULL:
cat("Time (secs) per simulation FULL: \n")
as.numeric(t4-t3,"secs")/N
# Times ratio
cat("Time ratio: \n")
as.numeric(t4-t3,"secs")/as.numeric(t2-t1,"secs")
#
save(list=ls(), file=paste(filename, "Rdata", sep="."))
sink()
#
