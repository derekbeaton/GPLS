rm(list=ls())

pls1_env <- new.env()
source('H:/Software/gpls/plscar_v1_simplified.R', pls1_env)

## within version tests
crossprod(pls1_env$LX)
t(pls1_env$LX) %*% pls1_env$LY
diag(t(pls1_env$LX) %*% pls1_env$LY) / pls1_env$Deltas
crossprod(pls1_env$Tmat)
cor(ca(pls1_env$Y.hat)$u, ca(pls1_env$Y.resid)$u)
cor(ca(pls1_env$Y.hat)$fi, ca(pls1_env$Y.resid)$fi)


pls2_env <- new.env()
source('H:/Software/gpls/plscar_v2_simplified.R', pls2_env)

## within version tests
crossprod(pls2_env$LX)
t(pls2_env$LX) %*% pls2_env$LY
diag(t(pls2_env$LX) %*% pls2_env$LY) / pls2_env$Deltas
crossprod(pls2_env$Tmat)
cor(ca(pls2_env$Y.hat)$u, ca(pls2_env$Y.resid)$u)
cor(ca(pls2_env$Y.hat)$fi, ca(pls2_env$Y.resid)$fi)


## between version tests
pls2_env$LX / pls1_env$LX
pls2_env$LY / pls1_env$LY
pls2_env$Deltas / pls1_env$Deltas
pls2_env$Betas / pls1_env$Betas
pls2_env$Tmat / pls1_env$Tmat
pls2_env$Y.hat / pls1_env$Y.hat
pls2_env$Y.resid / pls1_env$Y.resid


pls2_env$r2.x.cumulative / pls1_env$r2.x.cumulative
pls2_env$r2.y.cumulative / pls1_env$r2.y.cumulative


### time to level up to the next test.