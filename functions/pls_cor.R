pls_cor <- function(X, Y, center_x = TRUE, center_y = TRUE, scale_x = TRUE, scale_y = TRUE, components = 0, tol = .Machine$double.eps){

  
  X <- scale(X, center = center_x, scale = scale_x)
  Y <- scale(Y, center = center_y, scale = scale_y)
  
  gpls_cor(X, Y,  components = components, tol = tol)
  
}