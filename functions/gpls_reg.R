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
  
  stopped_early <- F
  
  if((components > X_rank) | (components < 1)){
    components <- X_rank
  }
  
  lx <- t_mat <- matrix(NA,nrow(X),components)
  ly <- matrix(NA,nrow(Y),components)
  
  predicted_u <- fi <- p <- u <- matrix(NA,ncol(X),components)  
  fj <- q <- v <- matrix(NA,ncol(Y),components)
  
  r2_x_cumulative <- r2_y_cumulative <- d <- betas <- rep(NA, components)
  
  X_reconstructeds <- X_residuals <- array(NA,dim=c(nrow(X),ncol(X),components))
  Y_reconstructeds <- Y_residuals <- array(NA,dim=c(nrow(Y),ncol(Y),components))
  
  X_deflate <- X
  Y_deflate <- Y
  
  for(i in 1:components){
    
    gplssvd_results <- gpls_cor(X_deflate, Y_deflate, 
                                XRW = XRW, YRW = YRW,
                                XLW = XLW, YLW = YLW,
                                components = 1, tol = tol)
    
    u[,i] <- gplssvd_results$u
    p[,i] <- gplssvd_results$p
    fi[,i] <- gplssvd_results$fi
    v[,i] <- gplssvd_results$v
    q[,i] <- gplssvd_results$q
    fj[,i] <- gplssvd_results$fj
    d[i] <- gplssvd_results$d
    lx[,i] <- gplssvd_results$lx
    ly[,i] <- gplssvd_results$ly
    
    t_mat[,i] <- lx[,i] / sqrt(sum(lx[,i]^2))
    betas[i] <- t(ly[,i]) %*% t_mat[,i]
    predicted_u[,i] <- t(t_mat[,i]) %*% ((XLW %^% (1/2)) %*% X_deflate %*% (XRW %^% (1/2)))
    
    
    X_reconstructeds[,,i] <- (XLW %^% (-1/2)) %*% (t_mat[,i] %o% predicted_u[,i]) %*% (XRW %^% (-1/2))
      X_reconstructeds[abs(X_reconstructeds) < tol] <- 0
    Y_reconstructeds[,,i] <- (YLW %^% (-1/2)) %*% ((t_mat[,i] * betas[i]) %o% v[,i]) %*% (YRW %^% (-1/2))
      Y_reconstructeds[abs(Y_reconstructeds) < tol] <- 0
    
    
    X_residuals[,,i] <- (X_deflate - X_reconstructeds[,,i])
      X_residuals[abs(X_residuals) < tol] <- 0
    Y_residuals[,,i] <- (Y_deflate - Y_reconstructeds[,,i])
      Y_residuals[abs(Y_residuals) < tol] <- 0
    
    X_deflate <- X_residuals[,,i]
    Y_deflate <- Y_residuals[,,i]
    
    
    r2_x_cumulative[i] <- (X_trace-sum( ( (XLW %^% (1/2)) %*%  X_deflate %*% (XRW %^% (1/2)) ) ^2)) / X_trace
    r2_y_cumulative[i] <- (Y_trace-sum( ( (YLW %^% (1/2)) %*%  Y_deflate %*% (YRW %^% (1/2)) ) ^2)) / Y_trace
    
    
    if( (sum(Y_deflate^2) < tol) & (i < components) ){
     
      stopped_early <- T
      warning("gpls_reg: Y is fully deflated. Stopping early.")
      
    }
    if( (sum(X_deflate^2) < tol) & (i < components) ){
      
      stopped_early <- T
      warning("gpls_reg: X is fully deflated. Stopping early.")
    
    }
    
    if(stopped_early){
      break
    }
    
  }
  
  if(stopped_early){
    u <- u[,1:i]
    p <- p[,1:i]
    fi <- fi[,1:i]
    v <- v[,1:i]
    q <- q[,1:i]
    fj <- fj[,1:i]
    d <- d[1:i]
    lx <- lx[,1:i]
    ly <- ly[,1:i]
    
    t_mat <- t_mat[,1:i]
    betas <- betas[1:i]
    predicted_u <- predicted_u[,1:i]
    
    X_reconstructeds <- X_reconstructeds[,,1:i]
    Y_reconstructeds <- Y_reconstructeds[,,1:i]
    X_residuals <- X_residuals[,,1:i]
    Y_residuals <- Y_residuals[,,1:i]
    
    r2_x_cumulative <- r2_x_cumulative[1:i]
    r2_y_cumulative <- r2_y_cumulative[1:i]
  }
  
  
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