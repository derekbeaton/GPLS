### need to check that weightedZX & weightedZY make for PLS-CA...

library(TExPosition)
library(ours)
library(GSVD)


### can we get back PLSC, CCA, and RRR from gplssvd? If so we're in some very serious business.

  ### from the GSVD help

# Three "two-table" technique examples
X <- scale(wine$objective)
Y <- scale(wine$subjective)
R <- t(X) %*% Y

### first GSVD vs GPLSSVD
### then against established alternate methods



#### PLS FIRST
## an example of partial least squares (correlation)
pls.res <- gsvd(R)
pls.res$lx <- X %*% pls.res$u
pls.res$ly <- Y %*% pls.res$v

gplssvd_pls <- GSVD:::gplssvd(X, Y)
teppls_res <- tepPLS(wine$objective, wine$subjective, graphs=F, scale1 = T, scale2 = T)


diag(t(pls.res$lx) %*% pls.res$ly)
diag(t(gplssvd_pls$lx) %*% gplssvd_pls$ly)
diag(t(teppls_res$TExPosition.Data$lx) %*% teppls_res$TExPosition.Data$ly)

gplssvd_pls$fj / teppls_res$TExPosition.Data$fj
gplssvd_pls$fi / teppls_res$TExPosition.Data$fi




cca.res <- gsvd(
  DAT=(crossprod(X) %^% -1) %*% R %*% (crossprod(Y) %^% -1),
  LW=crossprod(X),
  RW=crossprod(Y)
)

cca.res$maybe_lx <- (X %*% cca.res$u)
cca.res$maybe_ly <- (Y %*% cca.res$v)

cca.res$lx <- (X %*% cca.res$p)
cca.res$ly <- (Y %*% cca.res$q)



gplssvd_cca <- GSVD:::gplssvd(
  X, Y,
  # XLW=(tcrossprod(X)), 
  # YLW=(tcrossprod(Y)),
  XRW=crossprod(X),
  YRW=crossprod(Y)
)


cc.res <- cancor(X,Y)

(crossprod(X) %^% (1/2)) %*% (crossprod(X) %^% -1) %*% t(X) %*% Y %*% (crossprod(Y) %^% -1) %*% (crossprod(Y) %^% (1/2))



#### this needs to be solved:
(crossprod(X) %^% (1/2)) %*% (crossprod(X) %^% -1)



rrr.res <- gsvd(
  DAT=(crossprod(X) %^% -1) %*% R,
  LW=crossprod(X)
)
rrr.res$lx <- (X %*% rrr.res$p)
rrr.res$ly <- (Y %*% rrr.res$q)

(X %*% rrr.res$u)
(Y %*% rrr.res$v)






##### NOW PLSCA
data("snps.druguse")


plsca_res <- tepPLSCA(snps.druguse$DATA1, snps.druguse$DATA2, make_data1_nominal = T, make_data2_nominal = T, graphs=F)

DATA1_preproc <- ca.preproc(make.data.nominal(snps.druguse$DATA1))
DATA2_preproc <- ca.preproc(make.data.nominal(snps.druguse$DATA2))

this.R <- t(DATA1_preproc$weightedZx) %*% DATA2_preproc$weightedZx
svd.res <- svd(this.R)
LX <- DATA1_preproc$weightedZx %*% svd.res$u
LY <- DATA2_preproc$weightedZx %*% svd.res$v

### so this is the optimization.
diag(t(LX) %*% LY)[1:3]
### from here is where all the magic needs to happen.

## ok that worked well.
gplssvd_res <- gplssvd(DATA1_preproc$Zx, DATA2_preproc$Zx, 
                       1/DATA1_preproc$m, 1/DATA2_preproc$m,
                       1/DATA1_preproc$w, 1/DATA2_preproc$w)


