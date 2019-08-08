## simple plscar version 2: with the Zx & GPLSSVD matrices and, effectively, plain PLSR
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


ZX.in <- X_preproc$Zx
ZY.in <- Y_preproc$Zx

XLW <- diag(X_preproc$m)
  XLW_inv <- diag(1/X_preproc$m)
YLW <- diag(Y_preproc$m)
  YLW_inv <- diag(1/X_preproc$m)

XRW <- diag(X_preproc$w)
  XRW_inv <- diag(1/X_preproc$w)
YRW <- diag(Y_preproc$w)
  YRW_inv <- diag(1/Y_preproc$w)


## the loop
for(i in 1:comps){
  
  gplssvd_results <- gplssvd(X = ZX.in, XLW = XLW_inv, XRW = XRW_inv,
                             Y = ZY.in, YLW = YLW_inv, YRW = YRW_inv, 
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
  
  ### the algorithm itself.
  ### probably need more computational efficiency here...
  Tmat[,i] <- LX[,i] / sqrt(sum(LX[,i]^2))
  Betas[i] <- t(LY[,i]) %*% Tmat[,i]
    ### this is different because of, effectively, Ex or just the margins.
  pred.U[,i] <- t(Tmat[,i]) %*% ZX.in
  # pred.U[,i] <- t(Tmat[,i]) %*% sweep(sweep(ZX.in, 1, sqrt(X_preproc$m), "/"), 2, sqrt(X_preproc$w),"/")
  
  
  ###
  ## around here is where the issue lies. pred.U is by a constant, so X.recs is; and then so is Y.recs (unless) we switch to V[,i]
  ###
  X.recs[,,i] <- Tmat[,i] %o% pred.U[,i]
    X.recs[abs(X.recs) < tol] <- 0
  # X.hats[,,i] <- X.recs[,,i] * matrix(X.scale,nrow(X),ncol(X),byrow=T) + matrix(X.center,nrow(X),ncol(X),byrow=T)  
  # Y.recs[,,i] <- (Tmat[,i] * Betas[i]) %o% V[,i]#Q[,i]#(FJ[,i]/Deltas[i])
  # Y.recs[,,i] <- sweep(sweep(((Tmat[,i] * Betas[i]) %o% V[,i]),1,sqrt(Y_preproc$m),"*"),2,sqrt(Y_preproc$w),"*")
    Y.recs[,,i] <- (YLW %^% (1/2)) %*% ((Tmat[,i] * Betas[i]) %o% V[,i]) %*% (YRW %^% (1/2))
      ## these should use YLW_inv because that will be the general W...
      Y.recs[abs(Y.recs) < tol] <- 0
  # Y.hats[,,i] <- Y.recs[,,i] * matrix(Y.scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y.center,nrow(Y),ncol(Y),byrow=T)
  
  
    #####
    ### the problem needs to come back to here... it's about the residuals
    #####
  X.resids[,,i] <- (ZX.in - X.recs[,,i])
    X.resids[abs(X.resids) < tol] <- 0
  Y.resids[,,i] <- (ZY.in - Y.recs[,,i])
    Y.resids[abs(Y.resids) < tol] <- 0
  
  ZX.in <- X.resids[,,i]
  ZY.in <- Y.resids[,,i]
  
  ### the r2 shouldn't be cumulative if the effects aren't; it should be per component.
  ### the r2s here are not the same as in the other one; I need the weights.
  r2.x.cumulative[i] <- (X.trace-sum( ( (XLW_inv %^% (1/2)) %*%  ZX.in %*% (XRW_inv %^% (1/2)) ) ^2)) / X.trace
  r2.y.cumulative[i] <- (Y.trace-sum( ( (YLW_inv %^% (1/2)) %*%  ZY.in %*% (YRW_inv %^% (1/2)) ) ^2)) / Y.trace
    ## can these be re-done with ZX and ZY as is? I don't think so; requires SVD-able data matrix
  
  
  
}
# r2.x <- diff(c(0,r2.x.cumulative))
# r2.y <- diff(c(0,r2.y.cumulative))

# Y.rec <- sweep(sweep( Tmat %*% diag(Betas) %*% t(V),1,sqrt(Y_preproc$m),"*"),2,sqrt(Y_preproc$w),"*")
Y.rec <-  (YLW %^% (1/2)) %*% (Tmat %*% diag(Betas) %*% t(V)) %*% (YRW %^% (1/2))
Y.hat <- (Y.rec + Y_preproc$Ex) * sum(Y)
Y.resid <- ((Y_preproc$Zx - Y.rec) + Y_preproc$Ex) * sum(Y)



crossprod(LX)
t(LX) %*% LY
diag(t(LX) %*% LY) / Deltas
crossprod(Tmat)