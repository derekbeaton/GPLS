### actually this whole goddamn thing should be replaced with the gpls_reg as its core
  ## and then I have the bits and pieces around it
  ## so now go do the gpls_reg as the core... and then also gpls_cor as a core, too...

  ## the major drawback there is that I have to fiddle with those weights so things could slow down...
  

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
  
  ## preliminaries
  X_gsvd <- gsvd(X)
    X_rank <- length(X_gsvd$d)
    X_trace <- sum(X_gsvd$d^2)
  rm(X_gsvd)
  
  Y_gsvd <- gsvd(Y)
    Y_rank <- length(Y_gsvd$d)
    Y_trace <- sum(Y_gsvd$d^2)
  rm(Y_gsvd)
  
  
  if(components > X_rank){
    components <- X_rank
  }
  
  lx <- t_mat <- matrix(NA,nrow(X),components)
  ly <- matrix(NA,nrow(Y),components)
  
  predicted_u <- fi <- p <- u <- matrix(NA,ncol(X),components)  
  fj <- q <- v <- matrix(NA,ncol(Y),components)
  
  r2_x_cumulative <- r2_y_cumulative <- d <- betas <- rep(NA, components)
  
  X_hats <- X_reconstructeds <- X_residuals <- array(NA,dim=c(nrow(X),ncol(X),components_to_keep))
  Y_hats <- Y_reconstructeds <- Y_residuals <- array(NA,dim=c(nrow(Y),ncol(Y),components_to_keep))
  
  X_delfate <- X
  Y_deflate <- Y
  
  for(i in 1:components){
    
    gplssvd_results <- gplssvd(X_deflate, Y_deflate, k = 1)
    
    u[,i] <- gplssvd_results$u
    p[,i] <- gplssvd_results$p
    fi[,i] <- gplssvd_results$fi
    v[,i] <- gplssvd_results$v
    q[,i] <- gplssvd_results$q
    fj[,i] <- gplssvd_results$fj
    d[i] <- gplssvd_results$d
    lx[,i] <- gplssvd_results$lx
    ly[,i] <- gplssvd_results$ly
    
    
    t_mat[,i] <- LX[,i] / sqrt(sum(LX[,i]^2))
    betas[i] <- t(LY[,i]) %*% t_mat[,i]
    predicted_u[,i] <- t(t_mat[,i]) %*% X_deflate
    
    
    X_reconstructeds[,,i] <- t_mat[,i] %o% predicted_U[,i]
      X_reconstructed[abs(X_reconstructeds) < tol] <- 0
    Y_reconstructeds[,,i] <- (t_mat[,i] * betas[i]) %o% v[,i]
      Y_reconstructeds[abs(Y_reconstructeds) < tol] <- 0
    
    X_hats[,,i] <- X_reconstructeds[,,i] * matrix(X_scale,nrow(X),ncol(X),byrow=T) + matrix(X_center,nrow(X),ncol(X),byrow=T)
    Y_hats[,,i] <- Y_reconstructeds[,,i] * matrix(Y_scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y_center,nrow(Y),ncol(Y),byrow=T)

    X_residuals[,,i] <- (X_deflate - X_reconstructeds[,,i])
      X_residuals[abs(X_residuals) < tol] <- 0
    Y_residuals[,,i] <- (Y_deflate - Y_reconstructeds[,,i])
      Y_residuals[abs(Y_residuals) < tol] <- 0
    
    X_deflate <- X_residuals[,,i]
    Y_deflate <- Y_residuals[,,i]
    
    if(sum(X_delfate^2) < tol){
      r2_x_cumulative[i] <- 1
      break 
    }else{
      r2_x_cumulative[i] <- (X_trace - sum(X_delfate^2)) / X_trace
    }
    if(sum(Y_deflate^2)){
      r2_y_cumulative[i] <- 1
      break
    }else{
      r2_y_cumulative[i] <- (Y_trace - sum(Y_deflate^2)) / Y_trace
    }

  }
  
  Y_reconstructed <- t_mat %*% diag(betas) %*% t(v)
    Y_reconstructed[abs(Y_reconstructed) < tol] <- 0
  Y_residual <- Y - Y_reconstructed
  Y_hat <- Y_reconstructed * matrix(Y_scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y_center,nrow(Y),ncol(Y),byrow=T)
  
  
  return( list(
    d = d, u = u, v = v, lx = lx, ly = ly,
    p = p, q = q, fi = fi, fj = fj,
    t_mat = t_mat, predicted_u = predicted_u, betas = betas,
    X_reconstructeds = X_reconstructeds, X_residuals = X_residuals,
    Y_reconstructeds = Y_reconstructeds, Y_residuals = Y_residuals,
    X_hats = X_hats, Y_hats = Y_hats,
    r2_x = diff(c(0,r2_x_cumulative)), r2_y = diff(c(0,r2_y_cumulative)),
    Y_reconstructed = Y_reconstructed, Y_residual = Y_residual, Y_hat = Y_hat
  ) ) 
  
}