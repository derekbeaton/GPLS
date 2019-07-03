### need to check that weightedZX & weightedZY make for PLS-CA...

library(TExPosition)
library(ours)

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




### ok so what happens now?
source('./reframe/gplssvd.R')


  ## ok that worked well.
gplssvd_res <- gplssvd(DATA1_preproc$Zx, DATA2_preproc$Zx, 
                       1/DATA1_preproc$m, 1/DATA2_preproc$m,
                       1/DATA1_preproc$w, 1/DATA2_preproc$w)


