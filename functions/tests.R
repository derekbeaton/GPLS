## some tests to match against core/base techniques

rm(list=ls())

library(GSVD)
data("wine")
data("snps.druguse", package = "GSVD")
library(TExPosition)

source('./functions/utils.R')

source('./functions/gpls_cor.R')
source('./functions/rrr.R')
source('./functions/cca.R')
source('./functions/pls_cor.R')
source('./functions/plsca_cor.R')

source('./functions/gpls_reg.R')
source('./functions/pls_reg.R')
source('./functions/plsca_reg.R')

source('./functions/gpls_can.R')
source('./functions/pls_can.R')
source('./functions/plsca_can.R')


########################################
####################
##########
#####
##
# Some tests & conditions of all the ones in the "correlation" family (should be renamed...)
#           renaming b/c these are "single pass SVD" with no deflation
#   PLS-cor, PLSCA-cor, RRR/RDA, CCA
##
#####
##########
####################
########################################

### plsc
pls_cor_res <- pls_cor(wine$subjective, wine$objective)
texpo_plsc <- tepPLS(wine$subjective, wine$objective, scale1 = T, scale2 = T, graphs = F)

pls_cor_res$d / texpo_plsc$TExPosition.Data$pdq$Dv
pls_cor_res$lx / texpo_plsc$TExPosition.Data$lx
pls_cor_res$ly / texpo_plsc$TExPosition.Data$ly
pls_cor_res$fi / texpo_plsc$TExPosition.Data$fi
pls_cor_res$fj / texpo_plsc$TExPosition.Data$fj

### plsca
plsca_cor_res <- plsca_cor(make_data_disjunctive(snps.druguse$DATA1), make_data_disjunctive(snps.druguse$DATA2))
texpo_plsca <- tepPLSCA(snps.druguse$DATA1, snps.druguse$DATA2, make_data1_nominal = T, make_data2_nominal = T, graphs = F)

plsca_cor_res$d / texpo_plsca$TExPosition.Data$pdq$Dv
plsca_cor_res$lx / texpo_plsca$TExPosition.Data$lx
plsca_cor_res$ly / texpo_plsca$TExPosition.Data$ly
plsca_cor_res$fi / texpo_plsca$TExPosition.Data$fi
plsca_cor_res$fj / texpo_plsca$TExPosition.Data$fj


### cca
cca_res <- cca(wine$subjective, wine$objective, scale_x = F, scale_y = F)
cc.res <- cancor(wine$subjective, wine$objective)

cca_res$d / cc.res$cor
cc.res$xcoef[,1:5] / cca_res$p
cc.res$ycoef[,1:5] / cca_res$q
cca_res$d / diag(t(cca_res$lx) %*% cca_res$ly)

### see: http://ftp.uni-bayreuth.de/math/statlib/S/rrr.s
X <- scale(wine$subjective)
Y <- scale(wine$objective)
Sxy<-t(X)%*%Y
Sxmin12<-(crossprod(X) %^% (-1/2))
D<-Sxmin12%*%Sxy
svdobj<-svd(D)
U<-as.matrix(svdobj$u)
V<-as.matrix(svdobj$v)
L<-svdobj$d
Beta<-t(t(Sxmin12%*%U)*L)
Alpha<-V

rrr_res <- rrr(wine$subjective, wine$objective)
rda_res <- rrr(wine$subjective, wine$objective)


L / rrr_res$d
  # L / rda_res$d
U / rrr_res$u
V / rrr_res$v
Beta / rrr_res$beta_matrix

rrr_res$d / diag(t(rrr_res$lx) %*% rrr_res$ly)


paaaaause()

########################################
####################
##########
#####
##
# Some tests & conditions of the "regression" family (should be renamed)
#           renaming b/c these are multi-pass SVD with asymmetric deflation
#   PLS-reg, PLSCA-reg
##
#####
##########
####################
########################################
pls_reg_res <- pls_reg(wine$subjective, wine$objective, scale_x = T, scale_y = F)
gpls_reg_res <- gpls_reg(scale(wine$subjective, scale = T), scale(wine$objective, scale = F))


pls_reg_res$d / gpls_reg_res$d
pls_reg_res$lx / gpls_reg_res$lx
pls_reg_res$ly / gpls_reg_res$ly
pls_reg_res$fi / gpls_reg_res$fi
pls_reg_res$fj / gpls_reg_res$fj

pls_reg_res$d / diag(t(pls_reg_res$lx) %*% pls_reg_res$ly)
crossprod(pls_reg_res$u)
crossprod(pls_reg_res$tx)
crossprod(pls_reg_res$lx)

cor(pls_reg_res$Y_hat, pls_reg_res$Y_residual)
cor(tolerance.svd(pls_reg_res$Y_hat)$u, tolerance.svd(pls_reg_res$Y_residual)$u)


### now compare against ols
lm_res <- lm(as.matrix(wine$objective) ~ as.matrix(wine$subjective))

lm_res$fitted.values / pls_reg_res$Y_hat
lm_res$residuals / pls_reg_res$Y_residual

### now compare against pls::plsr
plsrplsr_res <- pls::plsr(as.matrix(wine$objective) ~ as.matrix(wine$subjective), ncomp = length(pls_reg_res$d), scale = T)

plsrplsr_res$fitted.values[,,length(pls_reg_res$d)] / pls_reg_res$Y_hat
plsrplsr_res$residuals[,,length(pls_reg_res$d)] / pls_reg_res$Y_residual

  ## what is wrong with these?
    #### nothing, it all depends on scale...
plsrplsr_res$scores / pls_reg_res$lx
  ## why are the LYs re-scaled by the d in pls::plsr? it sort of ruins the nice property of LX'LY
plsrplsr_res$Yscores / pls_reg_res$ly

colSums(pls_reg_res$u_hat * pls_reg_res$u_hat) / plsrplsr_res$Xvar

### now compare some of this against PLSC
pls_reg_res <- pls_reg(wine$subjective, wine$objective, scale_x = T, scale_y = T)
pls_reg_res$d[1] / pls_cor_res$d[1]
pls_reg_res$lx[,1] / pls_cor_res$lx[,1]
pls_reg_res$ly[,1] / pls_cor_res$ly[,1]
pls_reg_res$fi[,1] / pls_cor_res$fi[,1]
pls_reg_res$fj[,1] / pls_cor_res$fj[,1]


## plscar
  ### also needs an example like the final one in the preprint.
  ## needs some additional tests against lm()

plsca_reg_res <- plsca_reg(make_data_disjunctive(snps.druguse$DATA1), make_data_disjunctive(snps.druguse$DATA2))
## beginning of some tests:
# cor(tolerance.svd(ca.preproc(plsca_reg_res$Y_hat)$weightedZx)$u,tolerance.svd(ca.preproc(plsca_reg_res$Y_residual)$weightedZx)$u)

plsca_reg_res$d / diag(t(plsca_reg_res$lx) %*% plsca_reg_res$ly)
crossprod(plsca_reg_res$u)
crossprod(plsca_reg_res$tx)
crossprod(plsca_reg_res$lx)

cor(ca_preproc(plsca_reg_res$Y_hat)$Z, ca_preproc(plsca_reg_res$Y_residual)$Z)
cor(tolerance.svd(plsca_reg_res$Y_hat)$u, tolerance.svd(plsca_reg_res$Y_residual)$u)

### compare to lm()
  ## there are multiple ways to do this but this is fine for now.
plscar_lm <- lm(ca_preproc(make_data_disjunctive(snps.druguse$DATA2))$Z ~ ca_preproc(make_data_disjunctive(snps.druguse$DATA1))$Z)


plscar_lm$fitted.values / (ca_preproc(plsca_reg_res$Y_hat)$Z)
plscar_lm$residuals / (ca_preproc(plsca_reg_res$Y_residual)$Z)

#### against PLSCA

plsca_cor_res$d[1] / plsca_reg_res$d[1]
plsca_cor_res$lx[,1] / plsca_reg_res$lx[,1]
plsca_cor_res$ly[,1] / plsca_reg_res$ly[,1]
plsca_cor_res$fi[,1] / plsca_reg_res$fi[,1]
plsca_cor_res$fj[,1] / plsca_reg_res$fj[,1]





########################################
####################
##########
#####
##
# Some tests & conditions of the "canonical" family (should be renamed)
#           renaming b/c these are multi-pass SVD with symmetric deflation
#   PLS-can, PLSCA-can
##
#####
##########
####################
########################################

pls_can_res <- pls_can(wine$subjective, wine$objective, scale_x = T, scale_y = T)

pls_can_res$d / diag(t(pls_can_res$lx) %*% pls_can_res$ly)
crossprod(pls_can_res$u)
crossprod(pls_can_res$tx)
crossprod(pls_can_res$lx)

crossprod(pls_can_res$v)
crossprod(pls_can_res$ty)
crossprod(pls_can_res$ly)

# cor(pls_can_res$Y_hat, pls_can_res$Y_residual)
# cor(tolerance.svd(pls_can_res$Y_hat)$u, tolerance.svd(pls_can_res$Y_residual)$u)
  ## Y is empty here...
cor(pls_can_res$X_hat, pls_can_res$X_residual)
cor(tolerance.svd(pls_can_res$X_hat)$u, tolerance.svd(pls_can_res$X_residual)$u)

pls_reg_res$d[1] / pls_can_res$d[1]
pls_reg_res$lx[,1] / pls_can_res$lx[,1]
pls_reg_res$ly[,1] / pls_can_res$ly[,1]
pls_reg_res$fi[,1] / pls_can_res$fi[,1]
pls_reg_res$fj[,1] / pls_can_res$fj[,1]

  ## also good.
# pls_can_res <- pls_can(wine$subjective, wine$objective, scale_x = T, scale_y = T, components = 2)
# cor(pls_can_res$Y_hat, pls_can_res$Y_residual)
# cor(tolerance.svd(pls_can_res$Y_hat)$u, tolerance.svd(pls_can_res$Y_residual)$u)
# cor(pls_can_res$X_hat, pls_can_res$X_residual)
# cor(tolerance.svd(pls_can_res$X_hat)$u, tolerance.svd(pls_can_res$X_residual)$u)


plsca_can_res <- plsca_can(make_data_disjunctive(snps.druguse$DATA1), make_data_disjunctive(snps.druguse$DATA2))

plsca_can_res$d / diag(t(plsca_can_res$lx) %*% plsca_can_res$ly)
crossprod(plsca_can_res$u)
crossprod(plsca_can_res$tx)
crossprod(plsca_can_res$lx)

crossprod(plsca_can_res$v)
crossprod(plsca_can_res$ty)
crossprod(plsca_can_res$ly)


plsca_cor_res$d[1] / plsca_can_res$d[1]
plsca_cor_res$lx[,1] / plsca_can_res$lx[,1]
plsca_cor_res$ly[,1] / plsca_can_res$ly[,1]
plsca_cor_res$fi[,1] / plsca_can_res$fi[,1]
plsca_cor_res$fj[,1] / plsca_can_res$fj[,1]

