#' @export
#'
#' @title Partial least squares regression (PLSREG)
#'
#' @description Computes partial least squares "regression" between two data matrices by way of generalized PLS correlation
#'
#' @param X Data matrix with \emph{I} rows
#' @param Y Data matrix with \emph{I} rows
#' @param center_X For the \code{X} matrix: A parameter to pass through to \code{center} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{X}.
#' @param center_Y For the \code{Y} matrix: A parameter to pass through to \code{center} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{Y}.
#' @param scale_X For the \code{X} matrix: A parameter to pass through to \code{scale} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{X}.
#' @param scale_Y For the \code{Y} matrix: A parameter to pass through to \code{scale} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{Y}.
#' @param components The number of components to return. If < 1 then the maximum components will be returned. Default = 0.
#' @param tol default is .Machine$double.eps. A parameter to pass through to \code{\link[GSVD]{gplssvd}}; eliminates effectively zero, negative, or imaginary singular values (thus components).
#'
#'
#' @return A list of outputs based on \code{\link{gpls_reg}}
#' \item{d.orig}{A vector containing the singular values of DAT above the tolerance threshold (based on eigenvalues).}
#' \item{l.orig}{A vector containing the eigen values of DAT above the tolerance threshold (\code{tol}).}
#' \item{tau}{A vector that contains the (original) explained variance per component (via eigenvalues: \code{$l.orig}.}
#' \item{d}{A vector of length \code{min(length(d.orig), k)} containing the retained singular values}
#' \item{l}{A vector of length \code{min(length(l.orig), k)} containing the retained eigen values}
#' \item{u}{Left (rows) singular vectors. Dimensions are \code{nrow(DAT)} by k.}
#' \item{p}{Left (rows) generalized singular vectors.}
#' \item{fi}{Left (rows) component scores.}
#' \item{lx}{Latent variable scores for rows of \code{X}}
#' \item{v}{Right (columns) singular vectors.}
#' \item{q}{Right (columns) generalized singular vectors.}
#' \item{fj}{Right (columns) component scores.}
#' \item{lx}{Latent variable scores for rows of \code{Y}}
#'
#' @seealso \code{\link{gpls_reg}} \code{\link{pls_cor}} \code{\link{gpls_cor}} \code{\link[GSVD]{gplssvd}},
#'
#' @references
#' Abdi H., Eslami, A., Guillemot, V., & Beaton., D. (2018). Canonical correlation analysis (CCA). In R. Alhajj and J. Rokne (Eds.), \emph{Encyclopedia of Social Networks and Mining (2nd Edition).} New York: Springer Verlag.
#' Krishnan, A., Williams, L.J., McIntosh, A.R., & Abdi, H. (2011). Partial Least Squares (PLS) methods for neuroimaging: A tutorial and review. \emph{NeuroImage}, \bold{56}, 455-475.
#'
#' @examples
#'
#'  data("wine", package = "GSVD")
#'  plsreg_results <- pls_reg(wine$objective, wine$subjective)
#'
#'  \dontrun{
#'      ## like in RRR/RDA and CCA, PLSC maximizes the latent variables:
#'      diag( t(plsreg_results$lx) %*% plsreg_results$ly ) ## same as plsreg_results$d
#'      plsreg_results$d
#'
#'      ## but this is an asymmetric relationship:
#'      crossprod(plsreg_results$lx) #orthogonal
#'      t(plsreg_results$lx) %*% plsreg_results$ly ## a triangular matrix
#'  }
#'
#' @keywords multivariate, diagonalization, partial least squares

pls_reg <- function(X, Y, center_x = TRUE, center_y = TRUE, scale_x = TRUE, scale_y = TRUE, components = 0, tol = .Machine$double.eps){


  X <- scale(X, center = center_x, scale = scale_x)
    if(center_x){
      X_center <- attributes(X)$`scaled:center`
    }else{
      X_center <- rep(0, ncol(X))
    }
    if(scale_x){
      X_scale <- attributes(X)$`scaled:scale`
    }else{
      X_scale <- rep(1, ncol(X))
    }

  Y <- scale(Y, center = center_y, scale = scale_y)
    if(center_y){
      Y_center <- attributes(Y)$`scaled:center`
    }else{
      Y_center <- rep(0, ncol(Y))
    }
    if(scale_y){
      Y_scale <- attributes(Y)$`scaled:scale`
    }else{
      Y_scale <- rep(1, ncol(Y))
    }



  gpls_reg_results <- gpls_reg(X = X, Y = Y, components = components, tol = tol)

  gpls_reg_results$Y_hat <- gpls_reg_results$Y_reconstructed * matrix(Y_scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y_center,nrow(Y),ncol(Y),byrow=T)

  gpls_reg_results$X_hats <- array(NA,dim=c(nrow(X), ncol(X), length(gpls_reg_results$d)))
  gpls_reg_results$Y_hats <- array(NA,dim=c(nrow(Y), ncol(Y), length(gpls_reg_results$d)))
  for(i in 1:length(gpls_reg_results$d)){
    gpls_reg_results$X_hats[,,i] <- gpls_reg_results$X_reconstructeds[,,i] * matrix(X_scale,nrow(X),ncol(X),byrow=T) + matrix(X_center,nrow(X),ncol(X),byrow=T)
    gpls_reg_results$Y_hats[,,i] <- gpls_reg_results$Y_reconstructeds[,,i] * matrix(Y_scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y_center,nrow(Y),ncol(Y),byrow=T)
  }

  return(gpls_reg_results)

}
