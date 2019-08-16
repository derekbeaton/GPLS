gpls_can <- function(X, Y,
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
  
  if((components > X_rank) | (components > Y_rank) | (components < 1)){
    components <- min(X_rank, Y_rank)
  }
  
  lx <- matrix(NA,nrow(X),components)
  ly <- matrix(NA,nrow(Y),components)
  
  cx <- fi <- p <- u <- matrix(NA,ncol(X),components)  
  cy <- fj <- q <- v <- matrix(NA,ncol(Y),components)
  
  r2_x_cumulative <- r2_y_cumulative <- d <- rep(NA, components)
  
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
   
      ## and come back to this...
    cx[,i] <- t(X_deflate) %*% lx[,i]
    cy[,i] <- t(Y_deflate) %*% ly[,i]
    
    
    
      ## come back to this...
    # PX <- cx %*% solve(t(cx) %*% cx) %*% t(cx)
    # PY <- cy %*% solve(t(cy) %*% cy) %*% t(cy)
      ## and come back to this...
    X_reconstructeds[,,i] <- X_deflate %*% (cx[,i] %*% solve(t(cx[,i]) %*% cx[,i]) %*% t(cx[,i]))
      X_reconstructeds[abs(X_reconstructeds) < tol] <- 0
    Y_reconstructeds[,,i] <- Y_deflate %*% (cy[,i] %*% solve(t(cy[,i]) %*% cy[,i]) %*% t(cy[,i]))
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
    
    cx <- cx[,1:i]
    cy <- cy[,1:i]
    
    X_reconstructeds <- X_reconstructeds[,,1:i]
    Y_reconstructeds <- Y_reconstructeds[,,1:i]
    X_residuals <- X_residuals[,,1:i]
    Y_residuals <- Y_residuals[,,1:i]
    
    r2_x_cumulative <- r2_x_cumulative[1:i]
    r2_y_cumulative <- r2_y_cumulative[1:i]
  }
  
  
  return( list(
    d = d, u = u, v = v, lx = lx, ly = ly,
    p = p, q = q, fi = fi, fj = fj,
    X_reconstructeds = X_reconstructeds, X_residuals = X_residuals,
    Y_reconstructeds = Y_reconstructeds, Y_residuals = Y_residuals,
    r2_x = diff(c(0,r2_x_cumulative)), r2_y = diff(c(0,r2_y_cumulative)),
    cx= cx, cy = cy
  ) ) 
  
  
}