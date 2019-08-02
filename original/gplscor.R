gplscor <- function(ZX, ZY,  XLW = rep(1, nrow(X)), YLW = rep(1, nrow(Y)), XRW = rep(1, nrow(X)), YRW = rep(1, nrow(Y)), k = 0, tol = .Machine$double.eps){
  
  
  gplssvd_results <- gplssvd(ZX, ZY, XLW = XLW, YLW = YLW, XRW = XRW, YRW = YRW, k = k, tol = tol)
  return(gplssvd_results)

}