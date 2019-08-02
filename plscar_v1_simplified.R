## simple plscar version 1: with the weightedZx matrices and, effectively, plain PLSR
  ## listed by steps and using gplssvd()
  ## this works (for categorical; which is effectively a special case)
  ## but this needs generalization and use of the gplssvd() sextuplet.

## set up and needed stuff
library(GSVD)
library(ours)
library(ExPosition)

data("TRAITS")
data("SNPS")

tol <- .Machine$double.eps

X <- make.data.nominal(TRAITS)
Y <- make.data.nominal(SNPS)
X_preproc <- ca.preproc(X)
Y_preproc <- ca.preproc(Y)


## preliminaries
X.svd <- tolerance.svd(X_preproc$weightedZx)
X.rank <- length(X.svd$d)
X.trace <- sum(X.svd$d^2)
rm(X.svd)

Y.svd <- tolerance.svd(Y_preproc$weightedZx)
Y.rank <- length(Y.svd$d)
Y.trace <- sum(Y.svd$d^2)
rm(Y.svd)

comps <- X.rank


Y.isEmpty <- X.isEmpty <- F

LX <- Tmat <- matrix(0,nrow(X),comps)
LY <- matrix(0,nrow(Y),comps)

pred.U <- FI <- P <- U <- matrix(0,ncol(X),comps)  
FJ <- Q <- V <- matrix(0,ncol(Y),comps)

r2.x.cumulative <- r2.y.cumulative <- Deltas <- beta.vec <- Betas <- vector("numeric",comps)

X.hats <- X.recs <- X.resids <- array(0,dim=c(nrow(X),ncol(X),comps))
Y.hats <- Y.recs <- Y.resids <- array(0,dim=c(nrow(Y),ncol(Y),comps))


ZX.in <- X_preproc$weightedZx
ZY.in <- Y_preproc$weightedZx
# X_preproc$Zx / sweep(sweep(X_preproc$weightedZx, 1, sqrt(X_preproc$m), "*"), 2, sqrt(X_preproc$w), "*")



## the loop
for(i in 1:comps){
  
  gplssvd_results <- gplssvd(ZX.in, ZY.in, k = 1)
  
  U[,i] <- gplssvd_results$u
  P[,i] <- gplssvd_results$p
  FI[,i] <- gplssvd_results$fi
  V[,i] <- gplssvd_results$v
  Q[,i] <- gplssvd_results$q
  FJ[,i] <- gplssvd_results$fj
  Deltas[i] <- gplssvd_results$d
  LX[,i] <- gplssvd_results$lx
  LY[,i] <- gplssvd_results$ly
  
  ### the algorithm itself.
  ### probably need more computational efficiency here...
  Tmat[,i] <- LX[,i] / sqrt(sum(LX[,i]^2))
  Betas[i] <- t(LY[,i]) %*% Tmat[,i]
  pred.U[,i] <- t(Tmat[,i]) %*% ZX.in
  
  
  X.recs[,,i] <- Tmat[,i] %o% pred.U[,i]
    X.recs[abs(X.recs) < tol] <- 0
  X.hats[,,i] <- (sweep(sweep(X.recs[,,i], 1, sqrt(X_preproc$m), "*"), 2, sqrt(X_preproc$w), "*") + X_preproc$Ex) * sum(X)
  Y.recs[,,i] <- (Tmat[,i] * Betas[i]) %o% V[,i]
    Y.recs[abs(Y.recs) < tol] <- 0
  Y.hats[,,i] <- (sweep(sweep(Y.recs[,,i], 1, sqrt(Y_preproc$m), "*"), 2, sqrt(Y_preproc$w), "*") + Y_preproc$Ex) * sum(Y)
  
  
  X.resids[,,i] <- (ZX.in - X.recs[,,i])
    X.resids[abs(X.resids) < tol] <- 0
  Y.resids[,,i] <- (ZY.in - Y.recs[,,i])
    Y.resids[abs(Y.resids) < tol] <- 0
  
  ZX.in <- X.resids[,,i]
  ZY.in <- Y.resids[,,i]
  
  ### the r2 shouldn't be cumulative if the effects aren't; it should be per component.
  r2.x.cumulative[i] <- (X.trace-sum(ZX.in^2)) / X.trace
  r2.y.cumulative[i] <- (Y.trace-sum(ZY.in^2)) / Y.trace
  
}

r2.x <- diff(c(0,r2.x.cumulative))
r2.y <- diff(c(0,r2.y.cumulative))

Y.rec <- Tmat %*% diag(Betas) %*% t(V)
Y.hat <- (sweep(sweep(Y.rec, 1, sqrt(Y_preproc$m), "*"), 2, sqrt(Y_preproc$w), "*") + Y_preproc$Ex) * sum(Y)
Y.resid <- (sweep(sweep(Y_preproc$weightedZx - Y.rec, 1, sqrt(Y_preproc$m), "*"), 2, sqrt(Y_preproc$w), "*") + Y_preproc$Ex) * sum(Y)

crossprod(LX)
t(LX) %*% LY
diag(t(LX) %*% LY) / Deltas
crossprod(Tmat)
# cor(Y.hat, Y.resid)

  ## ok this works.
cor(ca(Y.hat)$u, ca(Y.resid)$u)
cor(ca(Y.hat)$fi, ca(Y.resid)$fi)