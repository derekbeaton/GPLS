plsca_cor <- function(X, Y, make_x_disjunctive = TRUE, make_y_disjunctive = TRUE, components = 0, tol = .Machine$double.eps){
 
  if(make_x_disjunctive){
    # X <-  
  }
  if(make_y_disjunctive){
    # Y <- 
  }
  
  X_ca_preproc <- ca_preproc(X)
  Y_ca_preproc <- ca_preproc(Y)
  
  gplssvd(X_ca_preproc$Z, Y_ca_preproc$Z,
          XLW = 1/X_ca_preproc$m, XRW = 1/X_ca_preproc$w,
          YLW = 1/Y_ca_preproc$m, YRW = 1/Y_ca_preproc$w,
          k = components, tol = tol)
  
  
}