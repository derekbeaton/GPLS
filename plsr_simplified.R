## simple plsr
  ## listed by steps and using gplssvd()

## set up and needed stuff
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


## preliminaries
X.svd <- tolerance.svd(X)
X.rank <- length(X.svd$d)
X.trace <- sum(X.svd$d^2)
rm(X.svd)

Y.svd <- tolerance.svd(Y)
Y.rank <- length(Y.svd$d)
Y.trace <- sum(Y.svd$d^2)
rm(Y.svd)

comps <- X.rank
# comps <- 2

Y.isEmpty <- X.isEmpty <- F

LX <- Tmat <- matrix(0,nrow(X),comps)
LY <- matrix(0,nrow(Y),comps)

pred.U <- FI <- P <- U <- matrix(0,ncol(X),comps)  
FJ <- Q <- V <- matrix(0,ncol(Y),comps)

r2.x.cumulative <- r2.y.cumulative <- Deltas <- beta.vec <- Betas <- vector("numeric",comps)

X.hats <- X.recs <- X.resids <- array(0,dim=c(nrow(X),ncol(X),comps))
Y.hats <- Y.recs <- Y.resids <- array(0,dim=c(nrow(Y),ncol(Y),comps))


ZX.in <- X
ZY.in <- Y


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
  X.hats[,,i] <- X.recs[,,i] * matrix(X.scale,nrow(X),ncol(X),byrow=T) + matrix(X.center,nrow(X),ncol(X),byrow=T)  
  Y.recs[,,i] <- (Tmat[,i] * Betas[i]) %o% V[,i]
    Y.recs[abs(Y.recs) < tol] <- 0
  Y.hats[,,i] <- Y.recs[,,i] * matrix(Y.scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y.center,nrow(Y),ncol(Y),byrow=T)
    
    
  X.resids[,,i] <- (ZX.in - X.recs[,,i])
    X.resids[abs(X.resids) < tol] <- 0
  Y.resids[,,i] <- (ZY.in - Y.recs[,,i])
    Y.resids[abs(Y.resids) < tol] <- 0
  
  ZX.in <- X.resids[,,i]
  ZY.in <- Y.resids[,,i]
  
  ### the r2 shouldn't be cumulative if the effects aren't; it should be per component.
  r2.x.cumulative[i] <- (X.trace-sum(ZX.in^2)) / X.trace
  r2.y.cumulative[i] <- (Y.trace-sum(ZY.in^2)) / Y.trace
  ## compute hats
  
}

r2.x <- diff(c(0,r2.x.cumulative))
r2.y <- diff(c(0,r2.y.cumulative))

## full model
Y.rec <- Tmat %*% diag(Betas) %*% t(V)
  Y.rec[abs(Y.rec) < tol] <- 0
Y.hat <-  Y.rec * matrix(Y.scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y.center,nrow(Y),ncol(Y),byrow=T)
Y.resid <- Y - Y.rec



plsr_res <- pls::plsr(as.matrix(wine$subjective) ~ as.matrix(wine$objective), ncomp = comps)
# plsdr_res <- plsdepot::plsreg2(as.matrix(wine$objective), as.matrix(wine$subjective), comps = X.rank)



## find the commonalities between pls::plsr & my code
colSums(pred.U * pred.U) / plsr_res$Xvar
Y.resid / plsr_res$residuals[,,comps]
Y.hat / plsr_res$fitted.values[,,comps]
plsr_res$scores / LX
## not quite...; but at least it is constants though I don't really know why.
plsr_res$Yscores / LY

pause()

## key relationships:
crossprod(LX)
t(LX) %*% LY
diag(t(LX) %*% LY) / Deltas
crossprod(Tmat)
cor(Y.hat, Y.resid)
cor(Y.rec, Y.resid)
# cor(Y.rec, Y.hat)

cor(tolerance.svd(Y.rec)$u, tolerance.svd(Y.resid)$u)




