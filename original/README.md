Generalized Partial Least Squares
================

Introduction
============

The generalized in generalized partial least squares (GPLS) is for both generalization of various PLS and canonical methods (e.g., canonical correlation) and also generalized across different data types (e.g., categorical, continuous, ordinal, mixed).

This is the (likely temporary, or transitionary) home for GPLS.

The current code for GPLS is *very* preliminary. Updates will be made over the coming weeks. A rudimentary illustration is below with some of the requried tests/conditions for various PLS algorithms.

Examples
========

The first is a simple set of tests.

``` r
require(TExPosition)
```

    ## Loading required package: TExPosition

    ## Loading required package: prettyGraphs

    ## Loading required package: ExPosition

``` r
require(GSVD)
```

    ## Loading required package: GSVD

``` r
require(ours)
```

    ## Loading required package: ours

    ## 
    ## Attaching package: 'ours'

    ## The following object is masked from 'package:ExPosition':
    ## 
    ##     expo.scale

``` r
source('gpls.R')
source('gplscor.R')
source('gplsreg.R')
source('gplscan.R')
source('plsca.R')

X <- makeNominalData(SNPS)
Y <- makeNominalData(TRAITS)


texpo.plsca <- tepPLSCA(X,Y,graphs=F)

gplsca.res_cor <- plsca(X,Y,pls.type="cor")
gplsca.res_reg <- plsca(X,Y,pls.type="reg")
gplsca.res_can <- plsca(X,Y,pls.type="can")
```

This shows equivalence between the "old" version in `TExPosition` and the new version here.

``` r
all.equal(gplsca.res_cor$d, texpo.plsca$TExPosition.Data$pdq$Dv)
```

    ## [1] TRUE

``` r
all(abs((abs(gplsca.res_cor$fj) - abs(texpo.plsca$TExPosition.Data$fj))) < sqrt(.Machine$double.eps))
```

    ## [1] TRUE

``` r
all(abs((abs(gplsca.res_cor$fi) - abs(texpo.plsca$TExPosition.Data$fi))) < sqrt(.Machine$double.eps))
```

    ## [1] TRUE

``` r
all(abs((abs(gplsca.res_cor$LX) - abs(texpo.plsca$TExPosition.Data$lx))) < sqrt(.Machine$double.eps))
```

    ## [1] TRUE

``` r
all(abs((abs(gplsca.res_cor$LY) - abs(texpo.plsca$TExPosition.Data$ly))) < sqrt(.Machine$double.eps))
```

    ## [1] TRUE

This shows equivalence between the first components of GPLSCA-COR, GPLSCA-REG, and GPLSCA-CAN.

``` r
all.equal(gplsca.res_cor$d[1], gplsca.res_reg$d[1])
```

    ## [1] TRUE

``` r
all.equal(gplsca.res_cor$d[1], gplsca.res_can$d[1])
```

    ## [1] TRUE

``` r
all(abs(abs(gplsca.res_cor$fi[,1]) - abs(gplsca.res_reg$fi[,1])) < sqrt(.Machine$double.eps))
```

    ## [1] TRUE

``` r
all(abs(abs(gplsca.res_cor$fi[,1]) - abs(gplsca.res_can$fi[,1])) < sqrt(.Machine$double.eps))
```

    ## [1] TRUE

``` r
all(abs(abs(gplsca.res_cor$fj[,1]) - abs(gplsca.res_reg$fj[,1])) < sqrt(.Machine$double.eps))
```

    ## [1] TRUE

``` r
all(abs(abs(gplsca.res_cor$fj[,1]) - abs(gplsca.res_can$fj[,1])) < sqrt(.Machine$double.eps))
```

    ## [1] TRUE

``` r
all(abs(abs(gplsca.res_cor$LX[,1]) - abs(gplsca.res_reg$LX[,1])) < sqrt(.Machine$double.eps))
```

    ## [1] TRUE

``` r
all(abs(abs(gplsca.res_cor$LX[,1]) - abs(gplsca.res_can$LX[,1])) < sqrt(.Machine$double.eps))
```

    ## [1] TRUE

``` r
all(abs(abs(gplsca.res_cor$LY[,1]) - abs(gplsca.res_reg$LY[,1])) < sqrt(.Machine$double.eps))
```

    ## [1] TRUE

``` r
all(abs(abs(gplsca.res_cor$LY[,1]) - abs(gplsca.res_can$LY[,1])) < sqrt(.Machine$double.eps))
```

    ## [1] TRUE

Finally the optimizations (orthogonality constraints).

``` r
all.equal(diag(t(gplsca.res_cor$LX) %*% gplsca.res_cor$LY), gplsca.res_cor$d)
```

    ## [1] TRUE

``` r
is.diagonal.matrix((t(gplsca.res_cor$LX) %*% gplsca.res_cor$LY))
```

    ## [1] TRUE

``` r
all.equal(diag(t(gplsca.res_reg$LX) %*% gplsca.res_reg$LY), gplsca.res_reg$d)
```

    ## [1] TRUE

``` r
  ## precision issues.
is.diagonal.matrix(abs(crossprod(gplsca.res_reg$LX)) > .Machine$double.eps)
```

    ## [1] TRUE

``` r
is.diagonal.matrix(abs(crossprod(gplsca.res_reg$LY)) > .Machine$double.eps)
```

    ## [1] FALSE

``` r
all.equal(diag(t(gplsca.res_can$LX) %*% gplsca.res_can$LY), gplsca.res_can$d)
```

    ## [1] TRUE

``` r
  ## precision issues.
is.diagonal.matrix(abs(crossprod(gplsca.res_can$LX)) > .Machine$double.eps)
```

    ## [1] TRUE

``` r
is.diagonal.matrix(abs(crossprod(gplsca.res_can$LX)) > .Machine$double.eps)
```

    ## [1] TRUE

Et voila.
