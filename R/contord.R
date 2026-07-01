#' @importFrom stats qt
#' @importFrom mvtnorm pmvt
#' @export

contord <- function (marginal, Sigma, support = list(), Spearman = FALSE, df=Inf, integerdf=TRUE, prob=FALSE)
{
  if (!all(unlist(lapply(marginal, function(x) (sort(x)==x & min(x)>0 & max(x)<1))))) stop("Error in assigning marginal distributions!")
  if(!isSymmetric(Sigma) | min(eigen(Sigma)$values)<0 | !all(diag(Sigma)==1)) stop("Correlation matrix not valid!")
  len <- length(support)
  k   <- length(marginal)
  kj  <- numeric(k)
  Sigmaord <- diag(k) # initialize the implied correlation matrix
  for (i in 1:k) {
    kj[i] <- length(marginal[[i]]) + 1
    if (len == 0) {
      support[[i]] <- 1:kj[i]
    }
    if (Spearman) {
      s1 <- c(marginal[[i]], 1)
      s2 <- c(0,marginal[[i]])
      support[[i]] <- (s1+s2)/2
    }
  }
  L <- vector("list", k)
  for (i in 1:k) {
    L[[i]] <- qt(marginal[[i]], df = df)
    L[[i]] <- c(-Inf, L[[i]], +Inf)     # thresholds
  }
  for (q in 1:(k - 1)) {
    for (r in (q + 1):k) {
      pij <- matrix(0, kj[q], kj[r])
      for (i in 1:kj[q]) {
        for (j in 1:kj[r]) {
          low <- rep(-Inf, k)
          upp <- rep(Inf, k)
          low[q] <- L[[q]][i]
          low[r] <- L[[r]][j]
          upp[q] <- L[[q]][i + 1]
          upp[r] <- L[[r]][j + 1]
          if(is.infinite(df) | (df==as.integer(df) & integerdf==TRUE))
          {
            pij[i, j] <- pmvt(low, upp, rep(0, k), corr = Sigma, df = df)
          } else {
            pij[i, j] <- pmvt.alt(low, upp, Sigma, df) # much more computationally expensive, but also more reliable
          }
          low <- rep(-Inf, k)
          upp <- rep(Inf, k)
        }
      }
      my <- sum(apply(pij, 2, sum) * support[[r]])
      sigmay <- sqrt(sum(apply(pij, 2, sum) * support[[r]]^2) - my^2)
      mx <- sum(apply(pij, 1, sum) * support[[q]])
      sigmax <- sqrt(sum(apply(pij, 1, sum) * support[[q]]^2) - mx^2)
      mij <- support[[q]] %*% t(support[[r]])
      muij <- sum(mij * pij)
      covxy <- muij - mx * my
      corxy <- covxy/(sigmax * sigmay) # implied (i,j) correlation
      Sigmaord[q, r] <- corxy
    }
  }
  SO <- as.matrix(forceSymmetric(Sigmaord))
  if(prob){ # construction of the k-variate probability table
    P <- array(NA, dim=kj)
    idx_grid <- do.call(expand.grid, c(lapply(kj, seq_len), KEEP.OUT.ATTRS = FALSE))
    colnames(idx_grid) <- paste0("j", seq_len(k))
    for (r in seq_len(nrow(idx_grid))) {
      idx <- as.integer(idx_grid[r, ])
      lower <- vapply(seq_len(k), function(j) L[[j]][ idx[j]     ], numeric(1))
      upper <- vapply(seq_len(k), function(j) L[[j]][ idx[j] + 1 ], numeric(1))
      # probability on k-dimensional hyper-rectangle
      if(is.infinite(df) | (df==as.integer(df) & integerdf==TRUE)) {
        p <- pmvt(lower = lower, upper = upper,
                  delta = rep(0,k), corr = Sigma,
                  df=df)
      } else {
        p <- pmvt.alt(lower, upper, Sigma, df=df)
      }
      P[matrix(idx, nrow = 1)] <- as.numeric(p)
    }
  }
  if(prob) return(list(pij=P,SigmaOrd=SO))
  else return(SO)
}
