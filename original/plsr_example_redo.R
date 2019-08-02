# plsr example; do-over with GSVD & GPLSSVD

## need some data
library(GSVD)
library(ours)
library(ExPosition)

data("wine")

tol <- .Machine$double.eps

X <- expo.scale(wine$objective, scale = F, center = T)
  X.center <- attributes(X)$`scaled:center`
  X.scale <- attributes(X)$`scaled:scale`
Y <- expo.scale(wine$subjective, scale = F, center = T)
  Y.center <- attributes(Y)$`scaled:center`
  Y.scale <- attributes(Y)$`scaled:scale`

## need some set up

X.svd <- tolerance.svd(X)
X.rank <- length(X.svd$d)
X.trace <- sum(X.svd$d^2)
rm(X.svd)

Y.svd <- tolerance.svd(Y)
Y.rank <- length(Y.svd$d)
Y.trace <- sum(Y.svd$d^2)
rm(Y.svd)

comps <- X.rank


LX <- Tmat <- matrix(0,nrow(X),comps)
LY <- matrix(0,nrow(Y),comps)

pred.U <- FI <- P <- U <- matrix(0,ncol(X),comps)  
FJ <- Q <- V <- matrix(0,ncol(Y),comps)

r2.x.cumulative <- r2.y.cumulative <- Deltas <- beta.vec <- Betas <- vector("numeric",comps)

X.hats <- X.recs <- X.resids <- array(0,dim=c(nrow(X),ncol(X),comps))
Y.hats <- Y.recs <- Y.resids <- array(0,dim=c(nrow(Y),ncol(Y),comps))


ZX.in <- X
ZY.in <- Y


for(i in 1:comps){
  
  gplssvd_results <- gplssvd(ZX.in, ZY.in,
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
  # Betas[i] <- t(Tmat[,i]) %*% (ZY.in %*% (FJ[,i] / Deltas[i]) )
  # beta.vec[i] <- t(LY[,i]) %*% Tmat[,i]
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
  
  ZX.in <- X.resids[,,i]
  ZY.in <- Y.resids[,,i]  
  
}

Y.rec <- Tmat %*% diag(Betas) %*% t(Q)
Y.resid <- Y - Y.rec

### need to pull up the HA version of PLSR that I translated.

plsr_res <- pls::plsr(as.matrix(wine$subjective) ~ as.matrix(wine$objective))

plsdr_res <- plsdepot::plsreg2(as.matrix(wine$objective), as.matrix(wine$subjective), comps = X.rank)
simpls_res <- plsdepot::simpls(as.matrix(wine$objective), as.matrix(wine$subjective), comps = X.rank)
