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
#source('/Professional/Software/ExPosition-Family/ExPosition2/SlimPosition/staging_code/sp.plsr.R')
source('../Code/Code/plsr.R')

data("vehicles")
X <- as.matrix(vehicles[,1:12])
Y <- as.matrix(vehicles[,13:16])

  ## I want all components some how.
pls.plsr <- pls::plsr(Y ~ X, data = data.frame(X=X,Y=Y)) ## this will tell 
simpls.res = simpls(X, Y,comps=pls.plsr$ncomp) # would be fine.
pd.plsr = plsreg2(X, Y, crosval=FALSE,comps = pls.plsr$ncomp)
sp.plsr.res <- plsr(X,Y)
lm.res <- lm(Y ~ X)



### ok now test some equivalencies...
pls.plsr$residuals[,,12] / lm.res$residuals ### final resids from pls:plsr are the residuals of OLS.
pls.plsr$fitted.values[,,12] / lm.res$fitted.values
diag(t(pls.plsr$scores) %*% pls.plsr$Yscores) ##only the upper diag is orthogonal; but I don't know what this is.
cor(pls.plsr$residuals[,,12],pls.plsr$fitted.values[,,12])


pd.plsr$resid / lm.res$residuals
pd.plsr$y.pred / lm.res$fitted.values
t(pd.plsr$x.scores) %*% pd.plsr$y.scores ##only the upper diag is orthogonal; but I don't know what this is.
cor(pd.plsr$resid,pd.plsr$y.pred)


sp.plsr.res$Y.resid / lm.res$residuals
sp.plsr.res$Y.hat / lm.res$fitted.values

t(sp.plsr.res$lx) %*% sp.plsr.res$ly ##only the upper diag is orthogonal
  ## because lx is orthogonal to itself but ly is not orthogonal to itself
cor(sp.plsr.res$Y.hat,sp.plsr.res$Y.resid)
  ### PLSR is fucking fascinating.
#cor(sp.plsr.res$X.hat,sp.plsr.res$X.resid)




## one last test: orthogonality of SVD results
hat.svd <- svd(expo.scale(sp.plsr.res$Y.hat,scale="SS1"))
resid.svd <- svd(expo.scale(sp.plsr.res$Y.resid,scale="SS1"))

cor(hat.svd$u,resid.svd$u) #orthogonal. this has to be a guarantee of what is returned from PLS-CA-R.



## just a stupid test.
X <- expo.scale(X,scale="SS1")
Y <- expo.scale(Y,scale="SS1")
proj.mat <- X %*% solve(t(X) %*% X) %*% t(X)
y.hat <- proj.mat %*% Y

svd.X <- svd(X)
test.y.hat <- svd.X$u %*% t(svd.X$u) %*% Y

wut <- y.hat * matrix(attributes(Y)$`scaled:scale`,nrow(Y),ncol(Y),byrow=T) + matrix(attributes(Y)$`scaled:center`,nrow(Y),ncol(Y),byrow=T)
wut / lm.res$fitted.values

test.wut <- test.y.hat * matrix(attributes(Y)$`scaled:scale`,nrow(Y),ncol(Y),byrow=T) + matrix(attributes(Y)$`scaled:center`,nrow(Y),ncol(Y),byrow=T)
test.wut / lm.res$fitted.values

  ## ok so I guess I get it now: OLS is y projected onto the (sort of) MD matrix.