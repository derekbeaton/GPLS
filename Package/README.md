Generalized partial least squares R package
================

To install the GPLS package please use `devtools`. GPLS has a dependency
on the GSVD package (found [here](https://github.com/derekbeaton/gsvd)).
Install as follows:

``` r
require(devtools)
devtools::install_github("derekbeaton/gsvd")
devtools::install_github("derekbeaton/gpls",subdir = "Package")
```
