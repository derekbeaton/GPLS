plsca_can <- function(X, Y, components = 0, tol = .Machine$double.eps){
  
  X_ca_preproc <- ca_preproc(X)
  Y_ca_preproc <- ca_preproc(Y)
  
  
  gpls_can_results <- gpls_can(X_ca_preproc$Z, Y_ca_preproc$Z,
                               XLW = diag(1/X_ca_preproc$m), XRW = diag(1/X_ca_preproc$w),
                               YLW = diag(1/Y_ca_preproc$m), YRW = diag(1/Y_ca_preproc$w),
                               components = components, tol = tol)
  
  
  gpls_can_results$X_residual <- ((X_ca_preproc$Z - gpls_can_results$X_reconstructed) + X_ca_preproc$E) * sum(X)
  gpls_can_results$X_hat <- (gpls_can_results$X_reconstructed + X_ca_preproc$E) * sum(X)
  
  gpls_can_results$Y_residual <- ((Y_ca_preproc$Z - gpls_can_results$Y_reconstructed) + Y_ca_preproc$E) * sum(Y)
  gpls_can_results$Y_hat <- (gpls_can_results$Y_reconstructed + Y_ca_preproc$E) * sum(Y)
  
  gpls_can_results$X_hats <- array(NA,dim=c(nrow(X), ncol(X), length(gpls_can_results$d)))
  gpls_can_results$Y_hats <- array(NA,dim=c(nrow(Y), ncol(Y), length(gpls_can_results$d)))
  for(i in 1:length(gpls_can_results$d)){
    gpls_can_results$X_hats[,,i] <- (gpls_can_results$X_reconstructeds[,,i] + X_ca_preproc$E) * sum(X)
    gpls_can_results$Y_hats[,,i] <- (gpls_can_results$Y_reconstructeds[,,i] + Y_ca_preproc$E) * sum(Y)
  }
  
  return(gpls_can_results)
  
  
}