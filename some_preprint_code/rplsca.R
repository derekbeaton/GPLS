
rplsca <- function(X,Y,pls.type=c("cor"),lambdas=seq(0,5,1),ridge.type="fast",comps=0, tol=.Machine$double.eps){

  rpls.res <- vector("list",length(lambdas))
  names(rpls.res) <- paste0("lambda_",lambdas)
  
  X.res <- ca.preproc(X)
  Y.res <- ca.preproc(Y)
  
  ## in formulation make this equality clear.
  # ZX <- diag(rowSums(X)) %*% (diag(1/X.res$m) %*% X.res$Zx)
  # ZY <- diag(rowSums(Y)) %*% (diag(1/Y.res$m) %*% Y.res$Zx)
  X.sum <- sum(X)
  Y.sum <- sum(Y)
  ZX <- X.res$Zx * X.sum
  ZY <- Y.res$Zx * Y.sum
  i <- 1
  for(lambda in lambdas){
  
    if(ridge.type=="fast"){
      XM <- rowSums(X) + (lambda * 1)
      XW <- colSums(X) + (lambda * 1)
      YM <- rowSums(Y) + (lambda * 1)
      YW <- colSums(Y) + (lambda * 1)
      
      rpls.res[[i]] <- gpls( sweep(ZX,1,sqrt(XM),'/'), sweep(ZY,1,sqrt(YM),'/'), 1/XW, 1/YW, comps = comps, pls.type=pls.type, tol=tol)
      
      ## simple rescaling to get it back to the exact specifications of MCA
      rpls.res[[i]]$fi <- rpls.res[[i]]$fi * sqrt(X.sum)
      rpls.res[[i]]$fj <- rpls.res[[i]]$fj * sqrt(Y.sum)
      rpls.res[[i]]$p <- rpls.res[[i]]$p / sqrt(X.sum)
      rpls.res[[i]]$q <- rpls.res[[i]]$q / sqrt(Y.sum)
      rpls.res[[i]]$pred.U <- rpls.res[[i]]$pred.U / sqrt(X.sum)
      rpls.res[[i]]$Y.rec <- rpls.res[[i]]$Y.rec / sqrt(Y.sum)
      rpls.res[[i]]$Y.resid <- rpls.res[[i]]$Y.resid / sqrt(Y.sum)
      rpls.res[[i]]$Y.recs <- rpls.res[[i]]$Y.recs / sqrt(Y.sum)
      rpls.res[[i]]$Y.resids <- rpls.res[[i]]$Y.resids / sqrt(Y.sum)
      rpls.res[[i]]$X.recs <- rpls.res[[i]]$X.recs / sqrt(X.sum)
      rpls.res[[i]]$X.resids <- rpls.res[[i]]$X.resids / sqrt(X.sum)
      
      
    }else{
      XM <- diag(rowSums(X)) + (lambda * matrix.generalized.inverse( tcrossprod(ZX) ))
      XW <- diag(colSums(X)) + (lambda * ( t(ZX) %*% matrix.generalized.inverse(tcrossprod(ZX)) %*% ZX ) )
      YM <- diag(rowSums(Y)) + (lambda * matrix.generalized.inverse( tcrossprod(ZY) ))
      YW <- diag(colSums(Y)) + (lambda * ( t(ZY) %*% matrix.generalized.inverse(tcrossprod(ZY)) %*% ZY ) )  
      
      rpls.res[[i]] <- gpls(((XM %^% (-1/2) %*% ZX)  ), ( (YM %^% (-1/2) %*% ZY)  ), XW %^% (-1), YW %^% (-1), comps = comps, pls.type=pls.type)
      
      ## simple rescaling to get it back to the exact specifications of MCA
      rpls.res[[i]]$fi <- rpls.res[[i]]$fi * sqrt(X.sum)
      rpls.res[[i]]$fj <- rpls.res[[i]]$fj * sqrt(Y.sum)
      rpls.res[[i]]$p <- rpls.res[[i]]$p / sqrt(X.sum)
      rpls.res[[i]]$q <- rpls.res[[i]]$q / sqrt(Y.sum)
      rpls.res[[i]]$pred.U <- rpls.res[[i]]$pred.U / sqrt(X.sum)
      rpls.res[[i]]$Y.rec <- rpls.res[[i]]$Y.rec / sqrt(Y.sum)
      rpls.res[[i]]$Y.resid <- rpls.res[[i]]$Y.resid / sqrt(Y.sum)
      rpls.res[[i]]$Y.recs <- rpls.res[[i]]$Y.recs / sqrt(Y.sum)
      rpls.res[[i]]$Y.resids <- rpls.res[[i]]$Y.resids / sqrt(Y.sum)
      rpls.res[[i]]$X.recs <- rpls.res[[i]]$X.recs / sqrt(X.sum)
      rpls.res[[i]]$X.resids <- rpls.res[[i]]$X.resids / sqrt(X.sum)
      
    }
    
    i <- i + 1
  }
  
  return(rpls.res)  
}