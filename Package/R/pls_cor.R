#'
#' @export
#'
#' @title Partial least squares correlation
#'
#' @description
#'
#' @param X
#' @param Y
#' @param center_X
#' @param center_Y
#' @param scale_X
#' @param scale_Y
#' @param components
#' @param tol
#'
#'
#' @return a list directly from \code{link{gpls_cor}}, which comes from \code{\link[GSVD]{gplssvd}}
#' \item{}{}
#'
#' @seealso
#'
#' @examples
#'
#' @keywords

pls_cor <- function(X, Y, center_x = TRUE, center_y = TRUE, scale_x = TRUE, scale_y = TRUE, components = 0, tol = .Machine$double.eps){


  X <- scale(X, center = center_x, scale = scale_x)
  Y <- scale(Y, center = center_y, scale = scale_y)

  gpls_cor(X, Y,  components = components, tol = tol)

}
