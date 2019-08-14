cca <- function(X, Y, k = 0, tol = .Machine$double.eps){

  X <- scale(X, center = center_x, scale = scale_x)
  Y <- scale(Y, center = center_y, scale = scale_y)
  
  gplssvd(
    X = X %^% (-1),
    Y = Y %^% (-1),
    XRW=crossprod(X),
    YRW=crossprod(Y)
  )
    
}