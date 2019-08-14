plsca_cor <- function(X, Y, components = 0, tol = .Machine$double.eps){
 
  X_ca_preproc <- ca_preproc(X)
  Y_ca_preproc <- ca_preproc(Y)
  
  gpls_cor(X_ca_preproc$Z, Y_ca_preproc$Z,
          XLW = 1/X_ca_preproc$m, XRW = 1/X_ca_preproc$w,
          YLW = 1/Y_ca_preproc$m, YRW = 1/Y_ca_preproc$w,
          components = components, tol = tol)
  
  
}