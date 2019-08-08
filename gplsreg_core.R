## gplsreg core function
## taken from plscar v2, but we don't need the two versions of each weight (e.g., YLW and YLW_inv)

## can this be used to do both PLSR & PLSCAR?
  ## also where can I gain speed ups/tricks?
  ## this will not be as optimized as it should be.

gpls_reg <- function(X, Y, 
                     XLW = diag(1, nrow(X)), YLW = diag(1, nrow(Y)),
                     XRW = diag(1, ncol(X)), YRW = diag(1, ncol(Y)), 
                     components_to_keep = 0, tol = .Machine$double.eps){
  
  
  
  X_gsvd <- gsvd(X, XLW, XRW)
  X_rank <- length(X_gsvd$d)
  X_trace <- sum(X_gsvd$d^2)
  rm(X_gsvd)
  
  Y_gsvd <- gsvd(Y, YLW, YRW)
  Y_rank <- length(Y_gsvd$d)
  Y_trace <- sum(Y_gsvd$d^2)
  rm(Y_gsvd)
  
  if(components_to_keep<=0){
    components_to_keep <- X_rank
  }
  
  # work these in later.
  # Y.isEmpty <- X.isEmpty <- F
  
  LX <- Tmat <- matrix(0,nrow(X),components_to_keep)
  LY <- matrix(0,nrow(Y),components_to_keep)
  
  FI <- P <- U <- pred_U <- matrix(0,ncol(X),components_to_keep)  
  FJ <- Q <- V <- matrix(0, ncol(Y), components_to_keep)
  
  r2.x.cumulative <- r2.y.cumulative <- Deltas <- Betas <- vector("numeric",components_to_keep)
  
  X_reconstructed <- X_residuals <- array(0,dim=c(nrow(X),ncol(X),components_to_keep))
  Y_reconstructed <- Y_residuals <- array(0,dim=c(nrow(Y),ncol(Y),components_to_keep))
  
  X_deflate <- X
  Y_deflate <- Y
  
  for(i in 1:components_to_keep){
    
    gplssvd_results <- gplssvd(X_deflate, Y_deflate, 
                             XLW = XLW, YLW = YLW, 
                             XRW = XRW, YRW = YRW,
                             k = 1)
    
    U[,i] <- gplssvd_results$u
    P[,i] <- gplssvd_results$p
    FI[,i] <- gplssvd_results$fi
    V[,i] <- gplssvd_results$v
    Q[,i] <- gplssvd_results$q
    FJ[,i] <- gplssvd_results$fj
    Deltas[i] <- gplssvd_results$d
    LX[,i] <- gplssvd_results$lx
    LY[,i] <- gplssvd_results$ly
    
    
    Tmat[,i] <- LX[,i] / sqrt(sum(LX[,i]^2))
    Betas[i] <- t(LY[,i]) %*% Tmat[,i]
    pred_U[,i] <- t(Tmat[,i]) %*% ((XLW %^% (1/2)) %*% X_deflate %*% (XRW %^% (1/2)))
    
    ### I think maybe this also deserves the weights...
    # X_reconstructed[,,i] <- Tmat[,i] %o% pred_U[,i]
    X_reconstructed[,,i] <- (XLW %^% (-1/2)) %*% (Tmat[,i] %o% pred_U[,i]) %*% (XRW %^% (-1/2))
      X_reconstructed[abs(X_reconstructed) < tol] <- 0
    Y_reconstructed[,,i] <- (YLW %^% (-1/2)) %*% ((Tmat[,i] * Betas[i]) %o% V[,i]) %*% (YRW %^% (-1/2))
      Y_reconstructed[abs(Y_reconstructed) < tol] <- 0
    
    # if the residuals are wrong, so is everything else...
    # X_residuals[,,i] <- (XLW %^% (-1/2)) %*% (X_deflate - X_reconstructed[,,i]) %*% (XRW %^% (-1/2))
    X_residuals[,,i] <- (X_deflate - X_reconstructed[,,i])
      X_residuals[abs(X_residuals) < tol] <- 0
    Y_residuals[,,i] <- (Y_deflate - Y_reconstructed[,,i])
      Y_residuals[abs(Y_residuals) < tol] <- 0
    
    X_deflate <- X_residuals[,,i]
    Y_deflate <- Y_residuals[,,i]
    
    
    r2.x.cumulative[i] <- (X_trace-sum( ( (XLW %^% (1/2)) %*%  X_deflate %*% (XRW %^% (1/2)) ) ^2)) / X_trace
    r2.y.cumulative[i] <- (Y_trace-sum( ( (YLW %^% (1/2)) %*%  Y_deflate %*% (YRW %^% (1/2)) ) ^2)) / Y_trace
    
  }
  
  ### I suppose this can be skipped and then this gpls_reg becomes the core of the functions.
  # Y.rec <-  (YLW %^% (-1/2)) %*% (Tmat %*% diag(Betas) %*% t(V)) %*% (YRW %^% (-1/2))
  ## here is a problem: not everything will be this way...
  # Y.hat <- (Y.rec + Y_preproc$Ex) * sum(Y)
  # Y.resid <- ((Y_preproc$Zx - Y.rec) + Y_preproc$Ex) * sum(Y)
  
  return( list(
    LX = LX, LY = LY, Deltas = Deltas,
    U = U, P = P, FI = FI,
    V = V, Q = Q, FJ = FJ,
    Tmat = Tmat, Betas = Betas, pred_U = pred_U,
    X_reconstructed = X_reconstructed, X_residuals = X_residuals,
    Y_reconstructed = Y_reconstructed, Y_residuals = Y_residuals,
    r2.x = diff(c(0,r2.x.cumulative)), r2.y = diff(c(0,r2.y.cumulative))
  ) ) 
  
}