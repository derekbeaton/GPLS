library(GSVD)
library(GPLS)
data("snps.druguse")

X <- make_data_disjunctive(snps.druguse$DATA1)
Y <- make_data_disjunctive(snps.druguse$DATA2)

ca_preproc_X <- ca_preproc(X)
  M_X <- diag(ca_preproc_X$m)
  W_X <- diag(ca_preproc_X$w)
  ZX_tilde <- (invsqrt_psd_matrix(M_X) %*% ca_preproc_X$Z %*% invsqrt_psd_matrix(W_X))  

ca_preproc_Y <- ca_preproc(Y)
  M_Y <- diag(ca_preproc_Y$m)
  W_Y <- diag(ca_preproc_Y$w)
  ZY_tilde <- (invsqrt_psd_matrix(M_Y) %*% ca_preproc_Y$Z %*% invsqrt_psd_matrix(W_Y))  

ZR <- t(ZX_tilde) %*% ZY_tilde  
plsca_reg_results <- plsca_reg(X, Y)




## illustration of maximization 

t(plsca_reg_results$lx) %*% plsca_reg_results$ly
diag(t(plsca_reg_results$lx) %*% plsca_reg_results$ly)
plsca_reg_results$d

t((invsqrt_psd_matrix(M_X)) %*% ca_preproc_X$Z %*% (t(MASS::ginv(W_X))) %*%  plsca_reg_results$p[,1]  ) %*%
  ((invsqrt_psd_matrix(M_Y)) %*% ca_preproc_Y$Z %*% (t(MASS::ginv(W_Y))) %*%  plsca_reg_results$q[,1]  )

t((invsqrt_psd_matrix(M_X)) %*% ca_preproc_X$Z %*% (invsqrt_psd_matrix(W_X)) %*% (invsqrt_psd_matrix(W_X)) %*%  plsca_reg_results$p[,1]  ) %*%
  ((invsqrt_psd_matrix(M_Y)) %*% ca_preproc_Y$Z %*% (invsqrt_psd_matrix(W_Y))%*% (invsqrt_psd_matrix(W_Y)) %*%  plsca_reg_results$q[,1]  )

t( ZX_tilde %*% (invsqrt_psd_matrix(W_X)) %*%  plsca_reg_results$p[,1]  ) %*%
  ( ZY_tilde %*% (invsqrt_psd_matrix(W_Y)) %*%  plsca_reg_results$q[,1]  )

t( ZX_tilde %*% (invsqrt_psd_matrix(W_X)) %*% (sqrt_psd_matrix(W_X)) %*%  plsca_reg_results$u[,1]  ) %*%
  ( ZY_tilde %*% (invsqrt_psd_matrix(W_Y)) %*% (sqrt_psd_matrix(W_Y)) %*%  plsca_reg_results$v[,1]  )

t(ZX_tilde %*% plsca_reg_results$u[,1]) %*% (ZY_tilde %*% plsca_reg_results$v[,1])

t(plsca_reg_results$u[,1]) %*% t(ZX_tilde) %*% ZY_tilde %*% plsca_reg_results$v[,1]

t(plsca_reg_results$u[,1]) %*% ZR %*% plsca_reg_results$v[,1]

## diagonal orthogonal.
t(plsca_reg_results$lx) %*% plsca_reg_results$lx
## diagonal orthogonal & identity
t(plsca_reg_results$p) %*% (t(MASS::ginv(W_X))) %*% plsca_reg_results$p
t(plsca_reg_results$tx) %*% plsca_reg_results$tx

## not orthogonal, but identity on the diagonal
t(plsca_reg_results$q) %*% (t(MASS::ginv(W_Y))) %*% plsca_reg_results$q



## illustration of double decomposition

(plsca_reg_results$tx %*% t(plsca_reg_results$u_hat)) / ZX_tilde
((invsqrt_psd_matrix(t(MASS::ginv(M_X)))) %*% (plsca_reg_results$tx %*% t(plsca_reg_results$u_hat)) %*% invsqrt_psd_matrix((t(MASS::ginv(W_X))))) / ca_preproc_X$Z


## sort of a redundant check as this is literally how it is defined.
Z_Y_hat <- invsqrt_psd_matrix(( t(MASS::ginv(M_Y)))) %*% (plsca_reg_results$tx %*% diag(plsca_reg_results$betas) %*% t(plsca_reg_results$v)) %*% invsqrt_psd_matrix((t(MASS::ginv(W_Y))))
Z_Y_hat / plsca_reg_results$Y_reconstructed
  ## also the same
Z_Y_hat / invsqrt_psd_matrix((t(MASS::ginv(M_Y)))) %*% ZX_tilde %*% ( t(MASS::ginv(plsca_reg_results$u_hat)) %*% diag(plsca_reg_results$betas) %*% t(plsca_reg_results$v) %*% invsqrt_psd_matrix((t(MASS::ginv(W_Y)))))


## illustration of OLS

## or if we were lazy, to obtain residualization:
# Z_Y_hat <- plsca_reg_results$Y_reconstructed

## these two are the same
ZX_tilde %*%  t(MASS::ginv( t(ZX_tilde) %*% ZX_tilde )) %*% t(ZX_tilde) %*% ZY_tilde
((plsca_reg_results$tx %*% diag(plsca_reg_results$betas) %*% t(plsca_reg_results$v)))

## these two are the same
ZX_tilde %*%  t(MASS::ginv(t(ZX_tilde) %*% ZX_tilde ))  %*% t(ZX_tilde) %*% ca_preproc_Y$Z
Z_Y_hat

