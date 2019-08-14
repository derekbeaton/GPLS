cca <- function(X, Y, center_x = TRUE, center_y = TRUE, scale_x = TRUE, scale_y = TRUE, components = 0, tol = .Machine$double.eps){
  
  
  X <- scale(X, center = center_x, scale = scale_x)
  Y <- scale(Y, center = center_y, scale = scale_y)
  
  ## a bit of a trick
  gpls_cor(
    X = X %^% (-1), Y = Y %^% (-1),
    XRW=crossprod(X), YRW=crossprod(Y), 
    components = components, tol = tol
  )
  
  
}