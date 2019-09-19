library(GSVD)
library(GPLS)
data("snps.druguse")

X <- make_data_disjunctive(snps.druguse$DATA1)
Y <- make_data_disjunctive(snps.druguse$DATA2)

ca_preproc_X <- ca_preproc(X)
  M_X <- diag(ca_preproc_X$m)
  W_X <- diag(ca_preproc_X$w)
  ZX_tilde <- (M_X %^% (-1/2)) %*% ca_preproc_X$Z %*% (W_X %^% (-1/2))  

ca_preproc_Y <- ca_preproc(Y)
  M_Y <- diag(ca_preproc_Y$m)
  W_Y <- diag(ca_preproc_Y$w)
  ZY_tilde <- (M_Y %^% (-1/2)) %*% ca_preproc_Y$Z %*% (W_Y %^% (-1/2))  

ZR <- t(ZX_tilde) %*% ZY_tilde  
plsca_reg_results <- plsca_reg(X, Y)




## illustration of maximization (Eq. 9)

t(plsca_reg_results$lx) %*% plsca_reg_results$ly
diag(t(plsca_reg_results$lx) %*% plsca_reg_results$ly)
plsca_reg_results$d

# t((M_X %^% (-1/2)) %*% ca_preproc_X$Z %*% (W_X %^% (-1)) %*% plsca_reg_results$p[,1]) %*%
#   ((M_Y %^% (-1/2)) %*% ca_preproc_Y$Z %*% (W_Y %^% (-1)) %*% plsca_reg_results$q[,1])

t((M_X %^% (-1/2)) %*% ca_preproc_X$Z %*% (W_X %^% (-1)) %*%  plsca_reg_results$p[,1]  ) %*%
  ((M_Y %^% (-1/2)) %*% ca_preproc_Y$Z %*% (W_Y %^% (-1)) %*%  plsca_reg_results$q[,1]  )


t((M_X %^% (-1/2)) %*% ca_preproc_X$Z %*% (W_X %^% (-1/2)) %*% (W_X %^% (1/2)) %*%  (W_X %^% (-1/2) %*% plsca_reg_results$u[,1])  ) %*%
  ((M_Y %^% (-1/2)) %*% ca_preproc_Y$Z %*% (W_Y %^% (-1/2))%*% (W_Y %^% (1/2)) %*%  (W_Y %^% (-1/2) %*% plsca_reg_results$v[,1])  )

t(ZX_tilde %*% plsca_reg_results$u[,1]) %*% (ZY_tilde %*% plsca_reg_results$v[,1])

t(plsca_reg_results$u[,1]) %*% t(ZX_tilde) %*% ZY_tilde %*% plsca_reg_results$v[,1]

t(plsca_reg_results$u[,1]) %*% ZR %*% plsca_reg_results$v[,1]

## diagonal identity, and orthogonal.
t(plsca_reg_results$p) %*% W_X %*% plsca_reg_results$p
t(plsca_reg_results$tx) %*% plsca_reg_results$tx
t(plsca_reg_results$lx) %*% plsca_reg_results$lx

## diagonal identity, but not orthogonal.
t(plsca_reg_results$q) %*% W_Y %*% plsca_reg_results$q



## illustration of double decomposition (Eq. 10)

(plsca_reg_results$tx %*% t(plsca_reg_results$u_hat)) / ZX_tilde

#
# (YLW %^% (-1/2)) %*% (tx %*% diag(betas) %*% t(v)) %*% (YRW %^% (-1/2))

((M_X %^% (-1)) %^% (-1/2)) %*% (plsca_reg_results$tx %*% diag(plsca_reg_results$betas) %*% t(plsca_reg_results$v)) %*% ((W_X %^% (-1)) %^% (-1/2)) / plsca_reg_results$Y_reconstructed
## and re-express it as




## illustration of OLS (Eq. 11)

ZX_tilde %*%  (( t(ZX_tilde) %*% ZX_tilde ) %^% (-1)) %*% t(ZX_tilde) %*% ca_preproc_Y$Z

ca_preproc_X$Z %*%  (( t(ca_preproc_X$Z) %*% ca_preproc_X$Z ) %^% (-1)) %*% t(ca_preproc_X$Z) %*% ca_preproc_X$Z

(ZX_tilde %*%  (( t(ZX_tilde) %*% ZX_tilde ) %^% (-1)) %*% t(ZX_tilde) %*% ca_preproc_Y$Z) / plsca_reg_results$Y_reconstructed
