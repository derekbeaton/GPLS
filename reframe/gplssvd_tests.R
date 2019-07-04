### need to check that weightedZX & weightedZY make for PLS-CA...

library(TExPosition)
library(ours)
library(GSVD)
library(rrr)

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


cca.res <- gsvd(
  DAT=R,
  LW=crossprod(X) %^% (-1),
  RW=crossprod(Y) %^% (-1)
)

cca.res$maybe_lx <- (X %*% cca.res$u)
cca.res$maybe_ly <- (Y %*% cca.res$v)

cca.res$lx <- (X %*% cca.res$p)
cca.res$ly <- (Y %*% cca.res$q)



gplssvd_cca <- GSVD:::gplssvd(
  X, Y,
  # XLW=(tcrossprod(X)), 
  # YLW=(tcrossprod(Y)),
  XRW=crossprod(X) %^% (-1),
  YRW=crossprod(Y) %^% (-1)
)


### this one is the answer.
gplssvd_cca <- GSVD:::gplssvd(
X = X %*% (crossprod(X) %^% (-1)),
Y = (Y %*% (crossprod(Y) %^% (-1))),
XRW=crossprod(X),
YRW=crossprod(Y)
)

gplssvd_cca2 <- GSVD:::gplssvd(
  X = X %^% (-1),
  Y = Y %^% (-1),
  XRW=crossprod(X),
  YRW=crossprod(Y)
)

  ### figure out all the optimizations here


cc.res <- cancor(X,Y)

svd.res0 <- svd((crossprod(X) %^% (1/2)) %*% (crossprod(X) %^% -1) %*% t(X) %*% Y %*% (crossprod(Y) %^% -1) %*% (crossprod(Y) %^% (1/2)))

svd.res <- svd((crossprod(X) %^% (-1/2)) %*% t(X) %*% Y %*% (crossprod(Y) %^% (-1/2)))


svd((crossprod(X) %^% (-1/2)) %*% t(X) %*% Y %*% (crossprod(Y) %^% (-1/2)))

(crossprod(X) %^% (1/2))
t(X %^% (1/2)) %*% (X %^% (1/2))


#### this needs to be solved:
(crossprod(X) %^% (1/2)) %*% (crossprod(X) %^% -1)


(crossprod(X) %^% (-1/2))


t(X %^% (1/2)) %*% (X %^% (-1))



### a good reference see rrr.nonmiss: http://ftp.uni-bayreuth.de/math/statlib/S/rrr.s

Sxy<-t(X)%*%Y
#Sxmin12<-sym.matr.power(t(X)%*%X,-0.5)
Sxmin12<-(crossprod(X) %^% (-1/2))
D<-Sxmin12%*%Sxy
svdobj<-svd(D)
U<-as.matrix(svdobj$u)
V<-as.matrix(svdobj$v)
L<-svdobj$d
Beta<-t(t(Sxmin12%*%U)*L)
Alpha<-V

rrr.res <- gsvd(
  DAT=(crossprod(X) %^% -1) %*% R,
  LW=crossprod(X)
)
rrr.res$lx <- (X %*% rrr.res$p)
rrr.res$ly <- (Y %*% rrr.res$q)

(X %*% rrr.res$u)
(Y %*% rrr.res$v)

### this seems to be the answer
gplssvd_rrr <- GSVD:::gplssvd(
  X = X %*% (crossprod(X) %^% (-1)),
  Y = Y,
  XRW=crossprod(X)
)

gplssvd_rrr2 <- GSVD:::gplssvd(
  X = X %^% (-1),
  Y = Y,
  XRW=crossprod(X)
)


rrr_res <- rrpack::rrr(Y,X)




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


