#'
#' @export
#'
#' @title Canonical correlation analysis
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

cca <- function(X, Y, center_X = TRUE, center_Y = TRUE, scale_X = TRUE, scale_Y = TRUE, components = 0, tol = .Machine$double.eps){


  X <- scale(X, center = center_X, scale = scale_X)
  Y <- scale(Y, center = center_Y, scale = scale_Y)

  ## a bit of a trick
  gpls_cor(
    X = X %^% (-1), Y = Y %^% (-1),
    XRW=crossprod(X), YRW=crossprod(Y),
    components = components, tol = tol
  )


}
