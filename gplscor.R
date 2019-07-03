gplscor <- function(ZX,ZY,WX=rep(1,ncol(ZX)),WY=rep(1,ncol(ZY)),comps=2, tol=.Machine$double.eps){
  
  gsvd.res <- gsvd(t(ZX) %*% ZY, LW = WX, RW = WY, k = comps, tol=tol)
  gsvd.res$LX <- ZX %*% (gsvd.res$fi %*% diag(1/gsvd.res$d))
  gsvd.res$LY <- ZY %*% (gsvd.res$fj %*% diag(1/gsvd.res$d))
  
  return(gsvd.res)
  
  
  
  
}