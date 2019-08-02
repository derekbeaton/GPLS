# plscar example; do-over with GSVD & GPLSSVD

## need some data & setup
library(GSVD)
library(ours)
library(ExPosition)

data("TRAITS")
data("SNPS")

tol <- .Machine$double.eps

X <- make.data.nominal(TRAITS)
Y <- make.data.nominal(SNPS)

## need some CA preprocessing

X_preproc <- ca.preproc(X)
Y_preproc <- ca.preproc(Y)

## need some set up

X.svd <- tolerance.svd(X_preproc$weightedZx)
X.rank <- length(X.svd$d)
X.trace <- sum(X.svd$d^2)
rm(X.svd)

Y.svd <- tolerance.svd(Y_preproc$weightedZx)
Y.rank <- length(Y.svd$d)
Y.trace <- sum(Y.svd$d^2)
rm(Y.svd)

comps <- X.rank

## some setup

LX <- Tmat <- matrix(0,nrow(X),comps)
LY <- matrix(0,nrow(Y),comps)

pred.U <- FI <- P <- U <- matrix(0,ncol(X),comps)  
FJ <- Q <- V <- matrix(0,ncol(Y),comps)

r2.x.cumulative <- r2.y.cumulative <- Deltas <- Betas <- vector("numeric",comps)

X.hats <- X.recs <- X.resids <- array(0,dim=c(nrow(X),ncol(X),comps))
Y.hats <- Y.recs <- Y.resids <- array(0,dim=c(nrow(Y),ncol(Y),comps))


ZX.in <- X_preproc$Zx
ZY.in <- Y_preproc$Zx


for(i in 1:comps){
  
  gplssvd_results <- gplssvd(ZX.in, ZY.in, 
                             XLW = 1/X_preproc$m, YLW = 1/Y_preproc$m, 
                             XRW = 1/X_preproc$w, YRW = 1/Y_preproc$w,
                             k = 1)
  
  U[,i] <- gplssvd_results$u
  P[,i] <- gplssvd_results$p
  FI[,i] <- gplssvd_results$fi
  V[,i] <- gplssvd_results$v
  Q[,i] <- gplssvd_results$q
  FJ[,i] <- gplssvd_results$fj
  Deltas[i] <- gplssvd_results$d
  LX[,i] <- gplssvd_results$lx
  LY[,i] <- gplssvd_results$ly
  
  ### probably need more computational efficiency here...
  Tmat[,i] <- LX[,i] / sqrt(sum(LX[,i]^2))
  Betas[i] <- t(Tmat[,i]) %*% (ZY.in %*% (FJ[,i] / Deltas[i]) )
  pred.U[,i] <- t(Tmat[,i]) %*% ZX.in
  
  
  ### probably need more computational efficiency here...
  Tmat[,i] <- LX[,i] / sqrt(sum(LX[,i]^2))
  Betas[i] <- t(LY[,i]) %*% Tmat[,i]
  pred.U[,i] <- t(Tmat[,i]) %*% ZX.in
  
  
  
  X.recs[,,i] <- Tmat[,i] %o% pred.U[,i]
    X.recs[abs(X.recs) < tol] <- 0
  Y.recs[,,i] <- (Tmat[,i] * Betas[i]) %o% Q[,i]
    Y.recs[abs(Y.recs) < tol] <- 0
  
  X.resids[,,i] <- (ZX.in - X.recs[,,i])
    X.resids[abs(X.resids) < tol] <- 0
  Y.resids[,,i] <- (ZY.in - Y.recs[,,i])
    Y.resids[abs(Y.resids) < tol] <- 0
  
    ### this won't be exactly right... 
  ZX.in <- X.resids[,,i]
  ZY.in <- Y.resids[,,i]  
  
    ### these aren't quite right because of above.
  ZX.svd.res <- tolerance.svd(ZX.in)
  r2.x.cumulative[i] <- (X.trace - sum(ZX.svd.res$d^2)) / X.trace
  ZY.svd.res <- tolerance.svd(ZY.in)
  r2.y.cumulative[i] <- (X.trace - sum(ZY.svd.res$d^2)) / Y.trace


}

Y.rec <- Tmat %*% diag(Betas) %*% t(Q)
Y.resid <- Y_preproc$Zx - Y.rec
  
  # return(
  #   list(
  #     fi = FI,
  #     fj = FJ,
  #     u = U,
  #     p = P,
  #     v = V,
  #     q = Q,
  #     d = Deltas,
  #     LX = LX,
  #     LY = LY,
  #     b = Betas,
  #     pred.U = pred.U,
  #     Tmat = Tmat,
  #     Y.resids = Y.resids,
  #     Y.resid = Y.resid,
  #     Y.recs = Y.recs,
  #     Y.rec = Y.rec,
  #     X.recs = X.recs,
  #     X.resids = X.resids,
  #     r2.x.cumulative=r2.x.cumulative,
  #     r2.x = diff(c(0,r2.x.cumulative)),
  #     r2.y.cumulative=r2.y.cumulative,
  #     r2.y = diff(c(0,r2.y.cumulative))
  #   )
  # )