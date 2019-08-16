## ugh.

rm(list=ls())

library(GSVD)
data("wine")

X <- scale(wine$subjective)
Y <- scale(wine$objective)


components <- nc <- min(dim(X), dim(Y))


Xold = X
Yold = Y
A = matrix(0, ncol(X), nc)
B = matrix(0, ncol(Y), nc)
Xt = matrix(0, nrow(X), nc)
Yu = matrix(0, nrow(X), nc)

Xresid <- Xrecon <- array(0, dim=c(nrow(X), ncol(X), nc))
Yresid <- Yrecon <- array(0, dim=c(nrow(Y), ncol(Y), nc))

for (h in 1:nc) {
  XY_svd = svd(t(Xold) %*% Yold)
  A[, h] = XY_svd$u[, 1]
  B[, h] = XY_svd$v[, 1]
  Xt[, h] = Xold %*% A[, h]
  Yu[, h] = Yold %*% B[, h]
  
  cx = t(Xold) %*% Xt[, h]
  Px = cx %*% solve(t(cx) %*% cx) %*% t(cx)
  cy = t(Yold) %*% Yu[, h]
  Py = cy %*% solve(t(cy) %*% cy) %*% t(cy)
  
  # Xold = Xold - (Xold %*% Px)
  # Yold = Yold - (Yold %*% Py)
  
  Xrecon[,,h] <- (Xold %*% Px)
  Xresid[,,h] <- (Xold - Xrecon[,,h])
  Xold <- Xresid[,,h]
  
  Yrecon[,,h] <- (Yold %*% Py)
  Yresid[,,h] <- (Yold - Yrecon[,,h])
  Yold <- Yresid[,,h]
}

  ## LX
x.scores = Xt
  ## U
x.wgs = A
  ## LY
y.scores = Yu
  ## V
y.wgs = B

simplsca_res <- simplsca(X, Y, comps = nc)


### this is actually a much much better approach than the standard simplsca
  ### I doubt this is a new approach, but, I'll need to look into this way

Xt_mat <- lx <- matrix(NA,nrow(X),components)
Yt_mat <- ly <- matrix(NA,nrow(Y),components)

Xpredicted_u <- cx <- fi <- p <- u <- matrix(NA,ncol(X),components)  
Ypredicted_u <- cy <- fj <- q <- v <- matrix(NA,ncol(Y),components)

d <- rep(NA, components)

X_reconstructeds2 <- X_residuals2 <- X_reconstructeds <- X_residuals <- array(NA,dim=c(nrow(X),ncol(X),components))
Y_reconstructeds2 <- Y_residuals2 <- Y_reconstructeds <- Y_residuals <- array(NA,dim=c(nrow(Y),ncol(Y),components))

X_deflate <- X
Y_deflate <- Y

tol <- .Machine$double.eps

for(i in 1:components){
  
  gplssvd_results <- gplssvd(X_deflate, Y_deflate, k = 1, tol = tol)
  
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
  # cx[,i] <- t(X_deflate) %*% lx[,i]
  # cy[,i] <- t(Y_deflate) %*% ly[,i]
  
  Xt_mat[,i] <- lx[,i] / sqrt(sum(lx[,i]^2))
  # betas[i] <- t(ly[,i]) %*% t_mat[,i]
  Xpredicted_u[,i] <- t(Xt_mat[,i]) %*% X_deflate
  
  Yt_mat[,i] <- ly[,i] / sqrt(sum(ly[,i]^2))
  # betas[i] <- t(ly[,i]) %*% t_mat[,i]
  Ypredicted_u[,i] <- t(Yt_mat[,i]) %*% Y_deflate
  
  ## come back to this...
  # PX <- cx %*% solve(t(cx) %*% cx) %*% t(cx)
  # PY <- cy %*% solve(t(cy) %*% cy) %*% t(cy)
  ## and come back to this...
  # X_reconstructeds[,,i] <- X_deflate %*% (cx[,i] %*% solve(t(cx[,i]) %*% cx[,i]) %*% t(cx[,i]))
  # Y_reconstructeds[,,i] <- Y_deflate %*% (cy[,i] %*% solve(t(cy[,i]) %*% cy[,i]) %*% t(cy[,i]))
  
  
  ### this is a far, far better solution to the "canonical" form of PLS.
    ## and easily adaptable under the generalized framework.
  
  X_reconstructeds[,,i] <- (Xt_mat[,i] %o% Xpredicted_u[,i])
  Y_reconstructeds[,,i] <- (Yt_mat[,i] %o% Ypredicted_u[,i])
  
  # X_residuals[,,i] <- (X_deflate - X_reconstructeds[,,i])
  # Y_residuals[,,i] <- (Y_deflate - Y_reconstructeds[,,i])
  
  X_residuals[,,i] <- (X_deflate - X_reconstructeds[,,i])
  Y_residuals[,,i] <- (Y_deflate - Y_reconstructeds[,,i])
  
  
  X_deflate <- X_residuals[,,i]
  Y_deflate <- Y_residuals[,,i]
  
}

source('functions/gpls_cor.R')
source('functions/gpls_reg.R')

gpls_reg_res <- gpls_reg(X, Y)
