

pls_reg <- function(X, Y, center_x = TRUE, center_y = TRUE, scale_x = TRUE, scale_y = TRUE, components = 0, tol = .Machine$double.eps){
  
  
  X <- scale(X, center = center_x, scale = scale_x)
    if(center_x){
      X_center <- attributes(X)$`scaled:center`
    }else{
      X_center <- rep(0, ncol(X))
    }
    if(scale_x){
      X_scale <- attributes(X)$`scaled:scale`
    }else{
      X_scale <- rep(1, ncol(X))
    }
  
  Y <- scale(Y, center = center_y, scale = scale_y)
    if(center_y){
      Y_center <- attributes(Y)$`scaled:center`
    }else{
      Y_center <- rep(0, ncol(Y))
    }
    if(scale_y){
      Y_scale <- attributes(Y)$`scaled:scale`
    }else{
      Y_scale <- rep(1, ncol(Y))
    }
  

  
  gpls_reg_results <- gpls_reg(X = X, Y = Y, components = components, tol = tol)
  
  
  gpls_reg_results$Y_reconstructed <- gpls_reg_results$t_mat %*% diag(gpls_reg_results$betas) %*% t(gpls_reg_results$v)
    gpls_reg_results$Y_reconstructed[abs(gpls_reg_results$Y_reconstructed) < tol] <- 0
  gpls_reg_results$Y_residual <- Y - gpls_reg_results$Y_reconstructed
  gpls_reg_results$Y_hat <- gpls_reg_results$Y_reconstructed * matrix(Y_scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y_center,nrow(Y),ncol(Y),byrow=T)
  
  
  gpls_reg_results$X_hats <- array(NA,dim=c(nrow(X), ncol(X), length(gpls_reg_results$d)))
  gpls_reg_results$Y_hats <- array(NA,dim=c(nrow(Y), ncol(Y), length(gpls_reg_results$d)))
  for(i in 1:length(gpls_reg_results$d)){
    gpls_reg_results$X_hats[,,i] <- gpls_reg_results$X_reconstructeds[,,i] * matrix(X_scale,nrow(X),ncol(X),byrow=T) + matrix(X_center,nrow(X),ncol(X),byrow=T)
    gpls_reg_results$Y_hats[,,i] <- gpls_reg_results$Y_reconstructeds[,,i] * matrix(Y_scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y_center,nrow(Y),ncol(Y),byrow=T)
  }
  
  return(gpls_reg_results)
  
}