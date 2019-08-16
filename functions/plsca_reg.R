plsca_reg <- function(X, Y, components = 0, tol = .Machine$double.eps){
  
  X_ca_preproc <- ca_preproc(X)
  Y_ca_preproc <- ca_preproc(Y)
  
  
  gpls_reg_results <- gpls_reg(X_ca_preproc$Z, Y_ca_preproc$Z,
                               XLW = diag(1/X_ca_preproc$m), XRW = diag(1/X_ca_preproc$w),
                               YLW = diag(1/Y_ca_preproc$m), YRW = diag(1/Y_ca_preproc$w),
                               components = components, tol = tol)
  
  gpls_reg_results$Y_reconstructed <- ( diag(1/Y_ca_preproc$m) %^% (-1/2)) %*% (gpls_reg_results$t_mat %*% diag(gpls_reg_results$betas) %*% t(gpls_reg_results$v)) %*% (diag(1/Y_ca_preproc$w) %^% (-1/2))
     gpls_reg_results$Y_reconstructed[abs(gpls_reg_results$Y_reconstructed) < tol] <- 0
  gpls_reg_results$Y_residual <- ((Y_ca_preproc$Z - gpls_reg_results$Y_reconstructed) + Y_ca_preproc$E) * sum(Y)
  gpls_reg_results$Y_hat <- (gpls_reg_results$Y_reconstructed + Y_ca_preproc$E) * sum(Y)
  
  # gpls_reg_results$X_hats <- array(NA,dim=c(nrow(X), ncol(X), length(gpls_reg_results$d)))
  # gpls_reg_results$Y_hats <- array(NA,dim=c(nrow(Y), ncol(Y), length(gpls_reg_results$d)))
  # for(i in 1:length(gpls_reg_results$d)){
  #   gpls_reg_results$X_hats[,,i] <- gpls_reg_results$X_reconstructeds[,,i] * matrix(X_scale,nrow(X),ncol(X),byrow=T) + matrix(X_center,nrow(X),ncol(X),byrow=T)
  #   gpls_reg_results$Y_hats[,,i] <- gpls_reg_results$Y_reconstructeds[,,i] * matrix(Y_scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y_center,nrow(Y),ncol(Y),byrow=T)
  # }
  
  return(gpls_reg_results)
  
  
}