### need to check that weightedZX & weightedZY make for PLS-CA...

library(TExPosition)
library(ours)
library(GSVD)
library(rrr)
library(rrpack)


### can we get back PLSC, CCA, and RRR from gplssvd? If so we're in some very serious business.

### from the GSVD help

# Three "two-table" technique examples
X <- scale(wine$subjective)
Y <- scale(wine$objective)


#### PLS
## an example of partial least squares (correlation)
gplssvd_pls <- GSVD:::gplssvd(X, Y)
teppls_res <- tepPLS(X, Y, graphs=F, center1 = F, center2 = F, scale1 = F, scale2 = F)


diag(t(gplssvd_pls$lx) %*% gplssvd_pls$ly)
diag(t(gplssvd_pls$lx) %*% gplssvd_pls$ly)
diag(t(teppls_res$TExPosition.Data$lx) %*% teppls_res$TExPosition.Data$ly)

gplssvd_pls$fj / teppls_res$TExPosition.Data$fj
gplssvd_pls$fi / teppls_res$TExPosition.Data$fi







#### CCA
cc.res <- cancor(X,Y)

gplssvd_cca <- GSVD:::gplssvd(
  X = X %*% (crossprod(X) %^% (-1)),
  Y = (Y %*% (crossprod(Y) %^% (-1))),
  XRW=crossprod(X),
  YRW=crossprod(Y)
)

cc.res$cor
gplssvd_cca$d.orig
diag(t(gplssvd_cca$lx) %*% gplssvd_cca$ly)

diag(t(cc.res$xcoef) %*% crossprod(X) %*% cc.res$xcoef)
diag(t(gplssvd_cca$p) %*% crossprod(X) %*% gplssvd_cca$p)

diag(t(cc.res$ycoef) %*% crossprod(Y) %*% cc.res$ycoef)
diag(t(gplssvd_cca$q) %*% crossprod(Y) %*% gplssvd_cca$q)


cc.res$xcoef
gplssvd_cca$p

cc.res$ycoef
gplssvd_cca$q

diag(t(gplssvd_cca$fi) %*% (crossprod(X) %^% (-1)) %*% gplssvd_cca$fi)
diag(t(gplssvd_cca$fj) %*% (crossprod(Y) %^% (-1)) %*% gplssvd_cca$fj)


## I don't think I know why this is the case...
gplssvd_cca2 <- GSVD:::gplssvd(
  X = X %^% (-1),
  Y = Y %^% (-1),
  XRW=crossprod(X),
  YRW=crossprod(Y)
)

  ## all equal; shows that this is a relatively easier way to compute CCA.
mapply(function(x,y){x/y}, gplssvd_cca, gplssvd_cca2)
  ### which turns into...


svd((crossprod(X) %^% (-1/2)) %*% (t(X) %*% Y) %*% (crossprod(Y) %^% (-1/2)))
  ## the immediate above is the Eq (1) from Cherry (1996). Then the "canonical vectors" are 
    # (crossprod(X) %^% (-1/2)) %*% u etc...
  ## one GSVD approach
svd((crossprod(X) %^% (1/2)) %*% (crossprod(X) %^% (-1)) %*% t(X) %*% Y %*% (crossprod(Y) %^% (1/2)) %*% (crossprod(Y) %^% (-1)))
  ## final, simplified GSVD approach.
svd((crossprod(X) %^% (1/2)) %*%  t(X %^% (-1)) %*% (Y %^% (-1)) %*% (crossprod(Y) %^% (1/2)))
    ### baller.
    ### we can see this through the cancellations of X or Y.
((crossprod(X) %^% (-1)) %*% t(X)) / (t(X %^% (-1)) %*% (X %^% (-1)) %*% t(X))
((crossprod(X) %^% (-1)) %*% t(X)) / t(X %^% (-1))



#### RRR / MLR / RDA
### a good reference see rrr.nonmiss: http://ftp.uni-bayreuth.de/math/statlib/S/rrr.s
Sxy<-t(X)%*%Y
Sxmin12<-(crossprod(X) %^% (-1/2))
D<-Sxmin12%*%Sxy
svdobj<-svd(D)
U<-as.matrix(svdobj$u)
V<-as.matrix(svdobj$v)
L<-svdobj$d
Beta<-t(t(Sxmin12%*%U)*L)
Alpha<-V


rrpack.rrr_res <- rrpack::rrr(Y,X)
rrr.rrr_res <- rrr::rrr(X,Y)



### this seems to be the answer
gplssvd_rrr <- GSVD:::gplssvd(
  X = X %*% (crossprod(X) %^% (-1)),
  Y = Y,
  XRW=crossprod(X)
)



L
gplssvd_rrr$d
diag(t(gplssvd_rrr$lx) %*% gplssvd_rrr$ly)

diag(t(gplssvd_rrr$p) %*% crossprod(X) %*% gplssvd_rrr$p)

diag(t(gplssvd_rrr$q) %*% gplssvd_rrr$q)


diag(t(gplssvd_rrr$fi) %*% (crossprod(X) %^% (-1)) %*% gplssvd_rrr$fi)
diag(t(gplssvd_rrr$fj)%*% gplssvd_rrr$fj)



gplssvd_rrr2 <- GSVD:::gplssvd(
  X = X %^% (-1),
  Y = Y,
  XRW=crossprod(X)
)


## all equal; shows that this is a relatively easier way to compute CCA.
mapply(function(x,y){x/y}, gplssvd_rrr, gplssvd_rrr2)
### which turns into...

(crossprod(X) %^% (1/2)) %*%  t(X %^% (-1)) %*% Y
svd((crossprod(X) %^% (1/2)) %*%  t(X %^% (-1)) %*% Y)
### baller.
### this is the same as above.





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
gplssvd_res <- GPLS::gplssvd(DATA1_preproc$Zx, DATA2_preproc$Zx, 
                       1/DATA1_preproc$m, 1/DATA2_preproc$m,
                       1/DATA1_preproc$w, 1/DATA2_preproc$w)



