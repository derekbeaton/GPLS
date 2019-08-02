plsca <- function(X, Y, pls.type=c("cor"), comps=0, tol=.Machine$double.eps){
  
  
  X.res <- ca.preproc(X)
  Y.res <- ca.preproc(Y)
 
  ZX.in <- sweep(X.res$Zx, 1, sqrt(X.res$m), "/")
  ZY.in <- sweep(Y.res$Zx, 1, sqrt(Y.res$m), "/")
  
  gplsrca.res <- gpls(ZX.in, ZY.in, 1/X.res$w, 1/Y.res$w, comps = comps, pls.type=pls.type, tol=tol)
  
  return(gplsrca.res)
   
}