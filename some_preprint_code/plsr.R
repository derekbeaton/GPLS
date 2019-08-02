plsr <- function(X, Y, center.X = T, scale.X = "SS1", center.Y = T, scale.Y = "SS1", k = 0, compact = T, tol=.Machine$double.eps){

  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows.")
  }

      ## maybe I need a replacement for these to use a sweep or apply
        ## and maybe I can updated expo.scale to be a more minimal "class"
  X.orig <- X <- expo.scale(X, scale = scale.X, center = center.X)
    X.center <- attributes(X)$`scaled:center`
    X.scale <- attributes(X)$`scaled:scale`
  Y.orig <- Y <- expo.scale(Y, scale = scale.Y, center = center.Y)
    Y.center <- attributes(Y)$`scaled:center`
    Y.scale <- attributes(Y)$`scaled:scale`

  X.svd <- tolerance.svd(X)
    X.trace <- sum(X.svd$d^2)
    X.rank <- length(X.svd$d)
    rm(X.svd)
  Y.svd <- tolerance.svd(Y)
    Y.trace <- sum(Y.svd$d^2)
    Y.rank <- length(Y.svd$d)
    rm(Y.svd)
  Y.isEmpty <- X.isEmpty <- F


  if(k<1 | k>X.rank){
    k <- X.rank #min(c(X.rank,Y.rank))
  }
  

  pred.u.mat <- fi.mat <- u.mat <- matrix(0,ncol(X),k)
  fj.mat <- v.mat <- matrix(0,ncol(Y),k)
  t.mat <- LX.mat <- LY.mat <- matrix(0,nrow(X),k)
  r2.x <- r2.y <- delta.vec <- beta.vec <- rep(0,k)
  if(!compact){
    X.resids <- X.hats <- array(0,dim=c(nrow(X),ncol(X),k))
  }
  Y.resids <- Y.hats <- array(0,dim=c(nrow(Y),ncol(Y),k))

  for(i in 1:k){
    print(i)
    
    res <- gsvd(t(X) %*% Y, k=1)
    u.mat[,i] <- res$u
    fi.mat[,i] <- res$fi
    v.mat[,i] <- res$v
    fj.mat[,i] <- res$fj
    delta.vec[i] <- res$d

    LX.mat[,i] <- X %*% as.matrix(u.mat[,i])
    LY.mat[,i] <- Y %*% as.matrix(v.mat[,i])

    ## what makes it PLSR, essentially.
    t.mat[,i] <- as.matrix(LX.mat[,i] / sqrt(sum(LX.mat[,i]^2)))
    pred.u.mat[,i] <- t(X) %*% t.mat[,i]
    beta.vec[i] <- c(t(LY.mat[,i]) %*% t.mat[,i])


    ## these are component-wise/iterative hats and resids...
    X.rec <- t.mat[,i] %*% t(pred.u.mat[,i])
    X <- X - X.rec
      X[which(abs(X) < tol)] <- 0
      r2.x[i] <- (X.trace-sum(X^2)) / X.trace
    if(!compact){
      X.hats[,,i] <- X.rec * matrix(X.scale,nrow(X),ncol(X),byrow=T) + matrix(X.center,nrow(X),ncol(X),byrow=T) ## can make more efficient through sweep
      X.resids[,,i] <- (X.orig * matrix(X.scale,nrow(X),ncol(X),byrow=T) + matrix(X.center,nrow(X),ncol(X),byrow=T)) - X.hats[,,i] ## can make more efficient through sweep
    }

    Y.rec <- t(beta.vec[i] * t(t.mat[,i] %*% t(v.mat[,i])))
    Y <- Y - Y.rec
      Y[which(abs(Y) < tol)] <- 0
      r2.y[i] <- (Y.trace-sum(Y^2)) / Y.trace
    Y.hats[,,i] <- Y.rec * matrix(Y.scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y.center,nrow(Y),ncol(Y),byrow=T) ## can make more efficient through sweep
    Y.resids[,,i] <- (Y.orig * matrix(Y.scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y.center,nrow(Y),ncol(Y),byrow=T)) - Y.hats[,,i] ## can make more efficient through sweep

    if(is.empty.matrix(X,tol=tol) | is.empty.matrix(Y,tol=tol) | is.empty.matrix(t(X) %*% Y,tol=tol)){
      break
    }
  }

    ## this is do-able for PLSC
  
    ## ok so this stops working when we don't go through all components...
  #Y.rec <- (t.mat * matrix(beta.vec,nrow(Y),i,byrow=T)) %*% t(v.mat)
  Y.rec <- (t.mat * matrix(beta.vec,nrow(Y),length(beta.vec),byrow=T)) %*% t(v.mat)
    Y.rec[which(abs(Y.rec) < tol)] <- 0
  Y.hat <-  Y.rec * matrix(Y.scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y.center,nrow(Y),ncol(Y),byrow=T)
    Y.hat[which(abs(Y.hat) < tol)] <- 0
  Y.resid <- (Y.orig * matrix(Y.scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y.center,nrow(Y),ncol(Y),byrow=T)) - Y.hat
  if(!compact){
    X.rec <- t.mat %*% t(pred.u.mat)
      X.rec[which(abs(X.rec) < tol)] <- 0
    X.hat <- X.rec * matrix(X.scale,nrow(X),ncol(X),byrow=T) + matrix(X.center,nrow(X),ncol(X),byrow=T)
      X.hat[which(abs(X.hat) < tol)] <- 0
    X.resid <- (X.orig * matrix(X.scale,nrow(X),ncol(X),byrow=T) + matrix(X.center,nrow(X),ncol(X),byrow=T)) - X.hat
  }

  if(compact){
    res <- list(fi=fi.mat, fj=fj.mat, d=delta.vec, u=u.mat, v=v.mat, lx=LX.mat, ly=LY.mat, Y.hats=Y.hats, Y.hat=Y.hat, Y.resid=Y.resid, r2.x=r2.x, r2.y=r2.y)
  }else{
    res <- list(fi=fi.mat, fj=fj.mat, d=delta.vec, u=u.mat, v=v.mat, lx=LX.mat, ly=LY.mat, beta=beta.vec, t=t.mat, pred.u=pred.u.mat, X.hats=X.hats, X.hat=X.hat, Y.hats=Y.hats, Y.hat=Y.hat, X.resids=X.resids, X.resid=X.resid, Y.resids=Y.resids, Y.resid=Y.resid, r2.x=r2.x, r2.y=r2.y)
  }

  return(res)

}
