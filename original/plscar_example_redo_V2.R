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

Tmat <- Tmat2 <- matrix(0,nrow(X),comps)
pred.U <- pred.U2 <- matrix(0,ncol(X),comps)  
r2.xs <- r2.ys <- r2.xs2 <- r2.ys2 <- Betas <- Betas2 <- vector("numeric",comps)

ZX.in <- X_preproc$Zx
ZY.in <- Y_preproc$Zx

another_ZX.in <- X_preproc$weightedZx
another_ZY.in <- Y_preproc$weightedZx

for(i in 1:comps){
  
  gplssvd_results <- gplssvd(ZX.in, ZY.in, 
                             XLW = 1/X_preproc$m, YLW = 1/Y_preproc$m, 
                             XRW = 1/X_preproc$w, YRW = 1/Y_preproc$w,
                             k = 1)
  
  ### probably need more computational efficiency here...
  Tmat[,i] <- gplssvd_results$lx / sqrt(sum(gplssvd_results$lx^2))
  # Betas[i] <- t(Tmat[,i]) %*% (ZY.in %*% (gplssvd_results$fj / gplssvd_results$d) )
  Betas[i] <- t(gplssvd_results$ly) %*% Tmat[,i]
  pred.U[,i] <- t(Tmat[,i]) %*% ZX.in
  
  
  ## I need to go back to plsr and simplify it greatly, then make use of the gplssvd() for it.
  # t(Tmat[,1]) %*% diag(sqrt(1/X_preproc$m)) %*% X_preproc$Zx %*% diag(sqrt(1/X_preproc$w)) / t(Tmat[,1]) %*% X_preproc$weightedZx
  
  ZX.in <- ZX.in - (Tmat[,i] %o% pred.U[,i])
  ZY.in <- ZY.in - ((Tmat[,i] * Betas[i]) %o% gplssvd_results$q[,1])
    
  
  
  gplssvd_results2 <- gplssvd(another_ZX.in, another_ZY.in,
                              k = 1)
  Tmat2[,i] <- gplssvd_results2$lx / sqrt(sum(gplssvd_results2$lx^2))
  Betas2[i] <- t(gplssvd_results2$ly) %*% Tmat2[,i]
  pred.U2[,i] <- t(Tmat2[,i]) %*% another_ZX.in
  
  another_ZX.in <- another_ZX.in - (Tmat2[,i] %o% pred.U2[,i])
  another_ZY.in <- another_ZY.in - ((Tmat2[,i] * Betas2[i]) %o% gplssvd_results$v[,1])
  
  if(all(abs(another_ZX.in) < tol)){
    r2.xs2[i] <- 1
  }else{
    ZX.svd.res <- tolerance.svd(another_ZX.in)
    r2.xs2[i] <- (X.trace - sum(ZX.svd.res$d^2)) / X.trace
  }
  
  ZY.svd.res <- tolerance.svd(another_ZY.in)
  r2.ys2[i] <- (Y.trace - sum(ZY.svd.res$d^2)) / Y.trace


}

# Y.rec <- Tmat %*% diag(Betas) %*% t(Q)
# Y.resid <- Y_preproc$Zx - Y.rec
#   
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

  ## this is my PLSR
test_plsr <- plsr(X_preproc$weightedZx, Y_preproc$weightedZx, F, F, F, F)

### in plain PLSR, isn't the original Y orthogonal to Y.resid?
  ### I really need to get back to this.
  ### step back a bit and do this with plain PLSR and via GPLSSVD()
cor(test_plsr$Y.hat, test_plsr$Y.resid)
cor(svd(test_plsr$Y.hat)$u, svd(test_plsr$Y.resid)$u)
cor(Y_preproc$weightedZx, test_plsr$Y.resid)

test_plsr$Y.resid / another_ZY.in
test_plsr$Y.resid - another_ZY.in




### I need a lot of clean up here...
  ### I really hvae a lot to re-organize here. I need to pull up the preprint version of plscar and plsr to make use of those too