

gpls_cor <- function(X, Y,
                  XLW = diag(nrow(X)), YLW = diag(nrow(Y)),
                  XRW = diag(ncol(X)), YRW = diag(ncol(Y)),
                  components = 0, tol = .Machine$double.eps){

    ## that's it.
    gplssvd(X, Y, XLW = XLW, YLW = YLW, XRW = XRW, YRW = YRW, k = components, tol = tol)

}
