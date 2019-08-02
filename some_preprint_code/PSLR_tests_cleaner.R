## plain PLSR test(s):
  ## run my PLSR (derivative of HA's)
  ## run GS's PLSR (SIMPLS or PLSREG2)
  ## run OLS
  ## run plsr from pls package

## What I want to test/check/understand:
  ## residuals of Y are identical across all
  ## Y is orthogonal to Yhat and/or Y resids.
  ## latent variable orthogonality properties...

rm(list=ls())

library(pls)
library(plsdepot)
library(GSVD)
library(ExPosition)
library(ours)
#source('/Professional/Software/ExPosition-Family/ExPosition2/SlimPosition/staging_code/sp.plsr.R')
#source('plsr.R')

data("vehicles")
X <- as.matrix(vehicles[,1:12])
Y <- as.matrix(vehicles[,13:16])

  ## I want all components some how.
pls.plsr <- pls::plsr(Y ~ X, data = data.frame(X=X,Y=Y)) ## this will tell 
pd.plsr = plsreg2(X, Y, crosval=FALSE,comps = pls.plsr$ncomp)
sp.plsr.res <- plsr(X,Y)
lm.res <- lm(Y ~ X)

## just a stupid test.
X2 <- expo.scale(X,scale="SS1")
Y2 <- expo.scale(Y,scale="SS1")
proj.mat <- X2 %*% solve(t(X2) %*% X2) %*% t(X2)
y.hat <- proj.mat %*% Y2
fitted.values <- y.hat * matrix(attributes(Y2)$`scaled:scale`,nrow(Y),ncol(Y),byrow=T) + matrix(attributes(Y2)$`scaled:center`,nrow(Y),ncol(Y),byrow=T)
residz <- Y - fitted.values


fitted.values / lm.res$fitted.values
sp.plsr.res$Y.hat / fitted.values
sp.plsr.res$Y.hat / lm.res$fitted.values
sp.plsr.res$Y.hat / pd.plsr$y.pred
sp.plsr.res$Y.hat / pls.plsr$fitted.values[,,12]


residz / lm.res$residuals
sp.plsr.res$Y.resid / residz
sp.plsr.res$Y.resid / lm.res$residuals
sp.plsr.res$Y.resid / pd.plsr$resid
sp.plsr.res$Y.resid / pls.plsr$residuals[,,12]

cor(sp.plsr.res$Y.hat, sp.plsr.res$Y.resid)

## one last test: orthogonality of SVD results
hat.svd <- svd(expo.scale(sp.plsr.res$Y.hat,scale="SS1"))
resid.svd <- svd(expo.scale(sp.plsr.res$Y.resid,scale="SS1"))
cor(hat.svd$u,resid.svd$u) #orthogonal. this has to be a guarantee of what is returned from PLS-CA-R.




# 
# svd.X <- svd(X2)
# test.y.hat <- svd.X$u %*% t(svd.X$u) %*% Y2
# test.y.hat2 <- svd.X$u %*% t(svd.X$u) %*% Y
# 
# test.wut <- test.y.hat * matrix(attributes(Y2)$`scaled:scale`,nrow(Y),ncol(Y),byrow=T) + matrix(attributes(Y2)$`scaled:center`,nrow(Y),ncol(Y),byrow=T)
# test.wut / lm.res$fitted.values
# 
#   ## ok so I guess I get it now: OLS is y projected onto the (sort of) MD matrix.


## ok now do the same but with PLSCAR
rm(list=ls())
source('plsr.R')
load(url("https://github.com/derekbeaton/PLSCA_Framework/blob/master/Data/SNPS.rda?raw=true"))
load(url("https://github.com/derekbeaton/PLSCA_Framework/blob/master/Data/TRAITS.rda?raw=true"))

X <- makeNominalData(SNPS)
Y <- makeNominalData(TRAITS)

X.res <- ca.preproc(X)
Y.res <- ca.preproc(Y)

  ## all the same...
proj.mat <- X.res$weightedZx %*% ((t(X.res$weightedZx) %*% X.res$weightedZx) %^% (-1)) %*% t(X.res$weightedZx)
# proj.mat2 <- X.res$weightedZx %*% matrix.generalized.inverse(t(X.res$weightedZx) %*% X.res$weightedZx) %*% t(X.res$weightedZx)
#  
# proj.mat3 <- X.res$Zx %*% ((t(X.res$Zx) %*% X.res$Zx) %^% (-1)) %*% t(X.res$Zx)
# proj.mat4 <- X.res$Zx %*% matrix.generalized.inverse(t(X.res$Zx) %*% X.res$Zx) %*% t(X.res$Zx)

  ## not the same...
y.hat <- proj.mat %*% Y.res$weightedZx
# y.hat2 <- proj.mat %*% Y.res$Zx
# fitted.values <- y.hat * matrix(attributes(Y2)$`scaled:scale`,nrow(Y),ncol(Y),byrow=T) + matrix(attributes(Y2)$`scaled:center`,nrow(Y),ncol(Y),byrow=T)
residz <- Y.res$weightedZx - y.hat

lm.res <- lm(Y.res$weightedZx ~ X.res$weightedZx)




plsr.res <- plsr(X.res$weightedZx, Y.res$weightedZx, center.X = F, scale.X = F, center.Y = F, scale.Y = F)
plsr.res2 <- plsr(X.res$Zx, Y.res$Zx, center.X = F, scale.X = F, center.Y = F, scale.Y = F)

y.hat / lm.res$fitted.values
plsr.res$Y.hat / y.hat
plsr.res$Y.hat / lm.res$fitted.values
plsr.res2$Y.hat / lm.res$fitted.values ## not the same...


residz / lm.res$residuals
plsr.res$Y.resid / residz
plsr.res$Y.resid / lm.res$residuals
plsr.res2$Y.resid / lm.res$residuals ## not the same...

  ## I forget exactly where these values come from...
1/(plsr.res2$Y.hat / lm.res$fitted.values)[1,]
1/(plsr.res2$Y.resid / lm.res$residuals)[1,]

### ok so next is to rebuild the data so that hat and resid stay orthogonal but are CA-friendly.
cor(plsr.res$Y.hat, plsr.res$Y.resid)

  ### these are "weighted" Zx. let's do Zx. t

## key thing here: I need to compute the SVDs of these two mtrices (hat & resid)
  ## then figure out the CAs from those
  ## start with the weightedZx ones, because that's actually closest to worst case
  ## then do so with the Zx
  ## and map between the two...
  ## then connect back to the gpls

## you had some other thoughts while running, but those are gone now...