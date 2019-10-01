library(GSVD)
library(GPLS)
data("snps.druguse")

X <- make_data_disjunctive(snps.druguse$DATA1)
Y <- make_data_disjunctive(snps.druguse$DATA2)

ca_preproc_X <- ca_preproc(X)
  M_X <- diag(ca_preproc_X$m)
  T_M_X <- diag(rowSums(X))
  W_X <- diag(ca_preproc_X$w)
  T_W_X <- diag(colSums(X))
  T_ZX <- ca_preproc_X$Z * sum(X)
  
ca_preproc_Y <- ca_preproc(Y)
  M_Y <- diag(ca_preproc_Y$m)
  T_M_Y <- diag(rowSums(Y))
  W_Y <- diag(ca_preproc_Y$w)
  T_W_Y <- diag(colSums(Y))
  T_ZY <- ca_preproc_Y$Z * sum(Y)


### plain 
plsca_cor_results <- plsca_cor(X, Y)

### Takane-ian (plain)
takane_plsca_cor_results <- gpls_cor(T_ZX, T_ZY, T_M_X  %^% (-1), T_M_Y  %^% (-1), T_W_X  %^% (-1), T_W_Y %^% (-1))

plsca_cor_results$d / takane_plsca_cor_results$d
plsca_cor_results$lx / takane_plsca_cor_results$lx
plsca_cor_results$ly / takane_plsca_cor_results$ly

  ## off by a constant scaling factor; which is just the sum(X) and sum(Y) respectively for the different sets
plsca_cor_results$fi / takane_plsca_cor_results$fi
plsca_cor_results$fj / takane_plsca_cor_results$fj

plsca_cor_results$p / takane_plsca_cor_results$p
plsca_cor_results$q / takane_plsca_cor_results$q


### Takane-ian (regularized at 0)
lambda <- 0

X_row_projection <- ((T_ZX %*% t(T_ZX)) %^% (-1)) * lambda
Y_row_projection <- ((T_ZY %*% t(T_ZY)) %^% (-1)) * lambda

X_col_projection <- (t(T_ZX) %*% X_row_projection %*% T_ZX) * lambda
Y_col_projection <- (t(T_ZY) %*% Y_row_projection %*% T_ZY) * lambda

takane_plsca_cor_results__lambda0 <- gpls_cor(T_ZX, T_ZY, 
                                              (T_M_X + X_row_projection)  %^% (-1), 
                                              (T_M_Y + Y_row_projection)  %^% (-1), 
                                              (T_W_X + X_col_projection) %^% (-1), 
                                              (T_W_Y + Y_col_projection) %^% (-1))

takane_plsca_cor_results$d / takane_plsca_cor_results__lambda0$d
takane_plsca_cor_results$lx / takane_plsca_cor_results__lambda0$lx
takane_plsca_cor_results$ly / takane_plsca_cor_results__lambda0$ly

takane_plsca_cor_results$fi / takane_plsca_cor_results__lambda0$fi
takane_plsca_cor_results$fj / takane_plsca_cor_results__lambda0$fj


### Takane-ian (regularized at 2)
lambda <- 2

X_row_projection <- ((T_ZX %*% t(T_ZX)) %^% (-1)) * lambda
Y_row_projection <- ((T_ZY %*% t(T_ZY)) %^% (-1)) * lambda

X_col_projection <- (t(T_ZX) %*% X_row_projection %*% T_ZX) * lambda
Y_col_projection <- (t(T_ZY) %*% Y_row_projection %*% T_ZY) * lambda

takane_plsca_cor_results__lambda2 <- gpls_cor(T_ZX, T_ZY, 
                                              (T_M_X + X_row_projection)  %^% (-1), 
                                              (T_M_Y + Y_row_projection)  %^% (-1), 
                                              (T_W_X + X_col_projection) %^% (-1), 
                                              (T_W_Y + Y_col_projection) %^% (-1))

takane_plsca_cor_results$d / takane_plsca_cor_results__lambda2$d
takane_plsca_cor_results$lx / takane_plsca_cor_results__lambda2$lx
takane_plsca_cor_results$ly / takane_plsca_cor_results__lambda2$ly

takane_plsca_cor_results$fi / takane_plsca_cor_results__lambda2$fi
takane_plsca_cor_results$fj / takane_plsca_cor_results__lambda2$fj



### Takane-ian regularized across lambdas with visualizations of fi & fj

X_row_projection <- ((T_ZX %*% t(T_ZX)) %^% (-1))
Y_row_projection <- ((T_ZY %*% t(T_ZY)) %^% (-1))

X_col_projection <- (t(T_ZX) %*% X_row_projection %*% T_ZX)
Y_col_projection <- (t(T_ZY) %*% Y_row_projection %*% T_ZY)

lambdas <- c(0,1,2,5,10,20)
lambda_cols <- viridisLite::cividis(length(lambdas))
counter <- 1

for(i in 1:length(lambdas)){

 
  
  takane_plsca_cor_results__reg <- gpls_cor(T_ZX, T_ZY, 
                                                (T_M_X + (X_row_projection * lambdas[i])) %^% (-1), 
                                                (T_M_Y + (Y_row_projection * lambdas[i])) %^% (-1), 
                                                (T_W_X + (X_col_projection * lambdas[i])) %^% (-1), 
                                                (T_W_Y + (Y_col_projection * lambdas[i])) %^% (-1))
  
  
  ### can change these out to $fj to also highlight the effect.
  if(i == 1){
    plot(takane_plsca_cor_results__reg$fi, col=lambda_cols[i], pch=20, axes = F, asp=NA, main="Takane regularization")
    abline(h=0)
    abline(v=0)
  }else{
    points(takane_plsca_cor_results__reg$fi, col=lambda_cols[i], pch=20)
    segments( last_scores[,1], last_scores[,2],  
              takane_plsca_cor_results__reg$fi[,1],takane_plsca_cor_results__reg$fi[,2],
              col = lambda_cols[i-1]
              )
  }

  last_scores <- takane_plsca_cor_results__reg$fi
  
}
legend("topright", legend = lambdas, col = lambda_cols, pch = 20, horiz = T)




### Hybrid Takane-ian/Allen-ian regularized approach 
Y_row_projection <- X_row_projection <- diag(nrow(X))

X_col_projection <- diag(ncol(X))
Y_col_projection <- diag(ncol(Y))

lambdas <- c(0,1,2,5,10,20)
lambda_cols <- viridisLite::cividis(length(lambdas))
counter <- 1

for(i in 1:length(lambdas)){
  
  
  
  takane_allen_plsca_cor_results__reg <- gpls_cor(T_ZX, T_ZY, 
                                            (T_M_X + (X_row_projection * lambdas[i])) %^% (-1), 
                                            (T_M_Y + (Y_row_projection * lambdas[i])) %^% (-1), 
                                            (T_W_X + (X_col_projection * lambdas[i])) %^% (-1), 
                                            (T_W_Y + (Y_col_projection * lambdas[i])) %^% (-1))
  
  
  ### can change these out to $fj to also highlight the effect.
  if(i == 1){
    plot(takane_allen_plsca_cor_results__reg$fi, col=lambda_cols[i], pch=20, axes = F, asp=NA, main="Takane/Allen 'truncated' regularization")
    abline(h=0)
    abline(v=0)
  }else{
    points(takane_allen_plsca_cor_results__reg$fi, col=lambda_cols[i], pch=20)
    segments( last_scores[,1], last_scores[,2],  
              takane_allen_plsca_cor_results__reg$fi[,1],takane_allen_plsca_cor_results__reg$fi[,2],
              col = lambda_cols[i-1]
    )
  }
  
  last_scores <- takane_allen_plsca_cor_results__reg$fi
  
}
legend("topright", legend = lambdas, col = lambda_cols, pch = 20, horiz = T)

