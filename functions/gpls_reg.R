gpls_reg <- function(X, Y,
                     XLW = diag(nrow(X)), YLW = diag(nrow(Y)), 
                     XRW = diag(ncol(X)), YRW = diag(ncol(Y)),
                     components = 0, tol = .Machine$double.eps){
  
  
  X_gsvd <- gsvd(X, XLW, XRW)
    X_rank <- length(X_gsvd$d)
    X_trace <- sum(X_gsvd$d^2)
  rm(X_gsvd)
  
  Y_gsvd <- gsvd(Y, YLW, YRW)
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
    
    # technically, we could (should) ship this call off to gpls_cor
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
    # predicted_u[,i] <- t(t_mat[,i]) %*% X_deflate
    predicted_u[,i] <- t(t_mat[,i]) %*% ((XLW %^% (1/2)) %*% X_deflate %*% (XRW %^% (1/2)))
    
    
    # X_reconstructeds[,,i] <- t_mat[,i] %o% predicted_U[,i]
    X_reconstructeds[,,i] <- (XLW %^% (-1/2)) %*% (t_mat[,i] %o% predicted_u[,i]) %*% (XRW %^% (-1/2))
      X_reconstructed[abs(X_reconstructeds) < tol] <- 0
    # Y_reconstructeds[,,i] <- (t_mat[,i] * betas[i]) %o% v[,i]
    Y_reconstructeds[,,i] <- (YLW %^% (-1/2)) %*% ((t_mat[,i] * betas[i]) %o% v[,i]) %*% (YRW %^% (-1/2))
      Y_reconstructeds[abs(Y_reconstructeds) < tol] <- 0
    
    
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
      # r2_x_cumulative[i] <- (X_trace - sum(X_delfate^2)) / X_trace
    }
    if(sum(Y_deflate^2)){
      r2_y_cumulative[i] <- 1
      break
    }else{
      # r2_y_cumulative[i] <- (Y_trace - sum(Y_deflate^2)) / Y_trace
    }
    
  }
  
  # Y_reconstructed <- t_mat %*% diag(betas) %*% t(v)
  Y_reconstructed <- (YLW %^% (-1/2)) %*% (t_mat %*% diag(betas) %*% t(v)) %*% (YRW %^% (-1/2))
    Y_reconstructed[abs(Y_reconstructed) < tol] <- 0
  Y_residual <- Y - Y_reconstructed
  
  
  return( list(
    d = d, u = u, v = v, lx = lx, ly = ly,
    p = p, q = q, fi = fi, fj = fj,
    t_mat = t_mat, predicted_u = predicted_u, betas = betas,
    X_reconstructeds = X_reconstructeds, X_residuals = X_residuals,
    Y_reconstructeds = Y_reconstructeds, Y_residuals = Y_residuals,
    r2_x = diff(c(0,r2_x_cumulative)), r2_y = diff(c(0,r2_y_cumulative)),
    Y_reconstructed = Y_reconstructed, Y_residual = Y_residual
  ) ) 
  
  
}