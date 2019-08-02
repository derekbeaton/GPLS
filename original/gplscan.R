gplscan <- function(ZX,ZY,WX=diag(ncol(ZX)),WY=diag(ncol(ZY)),comps=2, tol=.Machine$double.eps){
  
  break.out <- F
  
  X.svd <- gsvd(ZX,RW=WX,tol=tol)
  X.rank <- length(X.svd$d)
  X.trace <- sum(X.svd$d^2)
  rm(X.svd)
  Y.svd <- gsvd(ZY,RW=WY,tol=tol)
  Y.rank <- length(Y.svd$d)
  Y.trace <- sum(Y.svd$d^2)
  rm(Y.svd)
  
  if(comps > 0){
    comps <- min(c(comps,X.rank,Y.rank))
  }else{
    comps <- min(c(X.rank,Y.rank))
  }
  
  X.recs <- X.resids <- array(0,dim=c(nrow(ZX),ncol(ZX),comps))
  Y.recs <- Y.resids <- array(0,dim=c(nrow(ZY),ncol(ZY),comps))
  
  LX <- matrix(0,nrow(ZX),comps)
  LY <- matrix(0,nrow(ZY),comps)
  
  cx <- FI <- P <- U <- matrix(0,ncol(ZX),comps)  
  cy <- FJ <- Q <- V <- matrix(0,ncol(ZY),comps)
  
  r2.x.cumulative <- r2.y.cumulative <- Deltas <- vector("numeric",comps)
  
  ZX.in <- ZX
  ZY.in <- ZY
  
  for(i in 1:comps){
    
    gsvd.res <- try(gsvd(t(ZX.in) %*% ZY.in,WX,WY, k = 1, tol=tol))
    if(inherits(gsvd.res,"try-error")){
      #warning('GSVD::tolerance.svd() stopped. Results will be returned with fewer components than max(ncol(X),ncol(Y)).')
      break.out <- T
    }else{
      
      U[,i] <- gsvd.res$u
      P[,i] <- gsvd.res$p
      FI[,i] <- gsvd.res$fi
      V[,i] <- gsvd.res$v
      Q[,i] <- gsvd.res$q
      FJ[,i] <- gsvd.res$fj
      Deltas[i] <- gsvd.res$d
      LX[,i] <- ZX.in %*% as.matrix(gsvd.res$fi / gsvd.res$d) 
      LY[,i] <- ZY.in %*% as.matrix(gsvd.res$fj / gsvd.res$d)
      
      cx[,i] <- t(ZX.in) %*% LX[,i]
      cy[,i] <- t(ZY.in) %*% LY[,i]
      
      ## these are likely to blow up when things get big...
      ## can we come up with a P* that works like the regression version?
      # PX <- cx %*% solve(t(cx) %*% cx) %*% t(cx)
      # PY <- cy %*% solve(t(cy) %*% cy) %*% t(cy)
      
      # X.recs[,,i] <- (ZX.in %*% PX)
      X.recs[,,i] <- (ZX.in %*% (cx[,i] * (1/c(crossprod(cx[,i]))))) %*% cx[,i]
      X.recs[abs(X.recs) < tol] <- 0
      # Y.recs[,,i] <- (ZY.in %*% PY)
      Y.recs[,,i] <- (ZY.in %*% (cy[,i] * (1/c(crossprod(cy[,i]))))) %*% cy[,i]
      Y.recs[abs(Y.recs) < tol] <- 0
      
      X.resids[,,i] <- (ZX.in - X.recs[,,i])
      X.resids[abs(X.resids) < tol] <- 0
      Y.resids[,,i] <- (ZY.in - Y.recs[,,i])
      Y.resids[abs(Y.resids) < tol] <- 0
      
      ZX.in <- X.resids[,,i]
      ZY.in <- Y.resids[,,i]
      
      if( (sum(abs(ZX.in)) < tol) ){
        r2.x.cumulative[i] <- 1
        break.out <- T
      }else{
        ZX.gsvd.res <- gsvd(ZX.in,RW=WX,tol=tol)
        r2.x.cumulative[i] <- (X.trace - sum(ZX.gsvd.res$d^2)) / X.trace
        rm(ZX.gsvd.res)
      }
      
      if( (sum(abs(ZY.in)) < tol) ){
        r2.y.cumulative[i] <- 1
        break.out <- T
      }else{
        ZY.gsvd.res <- gsvd(ZY.in,RW=WY,tol=tol)
        r2.y.cumulative[i] <- (Y.trace - sum(ZY.gsvd.res$d^2)) / Y.trace
        rm(ZY.gsvd.res)
      }
      
      ## this may no longer be necessary... same as immediately above with the GSVD calls.
      test.r <- t(ZX.in) %*% ZY.in
      test.r[abs(test.r) < tol] <- 0
      if(sum(abs(test.r)) < tol){
        break.out <- T
      }
    }
    if(break.out){
      break
    }
    
  }
  
  return(
    list(
      fi = FI,
      fj = FJ,
      u = U,
      p = P,
      v = V,
      q = Q,
      d = Deltas,
      LX = LX,
      LY = LY,
      cx = cx,
      cy = cy,
      Y.resids = Y.resids,
      X.resids = X.resids,
      Y.recs = Y.recs,
      X.recs = X.recs,
      r2.x.cumulative=r2.x.cumulative,
      r2.x = diff(c(0,r2.x.cumulative)),
      r2.y.cumulative=r2.y.cumulative,
      r2.y = diff(c(0,r2.y.cumulative))
    )
  )  
}