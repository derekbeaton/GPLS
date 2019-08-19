#'
#' @export
#'
#' @title Reduced rank regression (see also redundancy analysis)
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


rrr <- function(X, Y, center_x = TRUE, center_y = TRUE, scale_x = TRUE, scale_y = TRUE, components = 0, tol = .Machine$double.eps){

  X <- scale(X, center = center_x, scale = scale_x)
  Y <- scale(Y, center = center_y, scale = scale_y)

  ## a bit of a trick
  rrr_res <- gpls_cor(
    X = X %^% (-1), Y = Y,
    XRW=crossprod(X),
    components = components, tol = tol
  )

  ## I should send back "Betas" here per http://ftp.uni-bayreuth.de/math/statlib/S/rrr.s
  rrr_res$beta_matrix <- (rrr_res$p %*% diag(rrr_res$d))

  rrr_res
}


#'
#' @export
#'
#' @title Redundancy analysis (see also reduced rank regression)
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

rda <- function(X, Y, center_x = TRUE, center_y = TRUE, scale_x = TRUE, scale_y = TRUE, components = 0, tol = .Machine$double.eps){

  rrr(X = X, Y = Y, center_x = center_x, center_y = center_y, scale_x = scale_x, scale_y = scale_y, components = components, tol = tol)

}
