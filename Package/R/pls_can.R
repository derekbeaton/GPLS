#' @export
#'
#' @title Partial least squares "canonical decomposition" (PLSCAN)
#'
#' @description Computes partial least squares "canonicaldecomposition " between two data matrices by way of generalized PLS canonical decomposition
#'
#' @param X Data matrix with \emph{I} rows and \emph{J} columns
#' @param Y Data matrix with \emph{I} rows and \emph{K} columns
#' @param center_X For the \code{X} matrix: A parameter to pass through to \code{center} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{X}.
#' @param center_Y For the \code{Y} matrix: A parameter to pass through to \code{center} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{Y}.
#' @param scale_X For the \code{X} matrix: A parameter to pass through to \code{scale} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{X}.
#' @param scale_Y For the \code{Y} matrix: A parameter to pass through to \code{scale} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{Y}.
#' @param components The number of components to return. If < 1 then the maximum components will be returned. Default = 0.
#' @param tol default is .Machine$double.eps. A parameter to pass through to \code{\link[GSVD]{gplssvd}}; eliminates effectively zero, negative, or imaginary singular values (thus components).
#'
#'
#' @return A list of outputs based on \code{\link{gpls_can}}
#' \item{d}{A vector containing the singular values from each iteration.}
#' \item{u}{Left (rows) singular vectors.}
#' \item{v}{Right (columns) singular vectors.}
#' \item{lx}{Latent variable scores for rows of \code{X}}
#' \item{ly}{Latent variable scores for rows of \code{Y}}
#' \item{p}{Left (rows) generalized singular vectors.}
#' \item{q}{Right (columns) generalized singular vectors.}
#' \item{fi}{Left (rows) component scores.}
#' \item{fj}{Right (columns) component scores.}
#' \item{tx}{\code{X} "Latent vectors": A normed version of \code{lx} for use in rebuilding \code{X} data}
#' \item{ty}{\code{Y} "Latent vectors": A normed version of \code{ly} for use in rebuilding \code{Y} data}
#' \item{u_hat}{\code{X} "Loading matrix": A "predicted" version of \code{u} for use in rebuilding \code{X} data}
#' \item{v_hat}{\code{Y} "Loading matrix": A "predicted" version of \code{v} for use in rebuilding \code{Y} data}
#' \item{X_reconstructeds}{A version of \code{X} reconstructed for each iteration (i.e., latent variable/component)}
#' \item{Y_reconstructeds}{A version of \code{Y} reconstructed for each iteration (i.e., latent variable/component)}
#' \item{X_residuals}{The residualized (i.e., \code{X - X_reconstructeds}) version of \code{X} for each iteration (i.e., latent variable/component)}
#' \item{Y_residuals}{The residualized (i.e., \code{Y - Y_reconstructeds}) version of \code{Y} for each iteration (i.e., latent variable/component)}
#' \item{r2_x}{Proporition of explained variance from \code{X} to each latent variable/component.}
#' \item{r2_y}{Proporition of explained variance from \code{Y} to each latent variable/component.}
#' \item{X_reconstructed}{A version of \code{X} reconstructed from all iterations (i.e., latent variables/components); see \code{components}.}
#' \item{X_residual}{The residualized (i.e., \code{X - X_reconstructed} from all iterations (i.e., latent variables/components); see \code{components}.}
#' \item{Y_reconstructed}{A version of \code{Y} reconstructed from all iterations (i.e., latent variables/components); see \code{components}.}
#' \item{Y_residual}{The residualized (i.e., \code{Y - Y_reconstructed} from all iterations (i.e., latent variables/components); see \code{components}.}
#' \item{X_hat}{The re-centered and re-scaled version of \code{X_reconstructed} to have the same center and scale as \code{X}.}
#' \item{Y_hat}{The re-centered and re-scaled version of \code{Y_reconstructed} to have the same center and scale as \code{Y}.}
#' \item{X_hats}{The re-centered and re-scaled versions of \code{X_reconstructeds} per iteration (i.e., latent variable/component) to have the same center and scale as \code{X}.}
#' \item{Y_hats}{The re-centered and re-scaled versions of \code{Y_reconstructeds} per iteration (i.e., latent variable/component) to have the same center and scale as \code{Y}.}
#'
#'
#' @seealso \code{\link{plsca_can}} \code{\link{gpls_can}} \code{\link[GSVD]{gplssvd}}
#'
#' @references
#' Tenenhaus, M. (1998). La Regression PLS. Theorie et Pratique. \emph{Editions TECHNIP}, Paris.
#'
#'
#' @examples
#'  data("wine", package = "GSVD")
#'  plscan_results <- pls_can(wine$objective, wine$subjective)
#'
#'  \dontrun{
#'      ## PLSCOR, PLSCAN maximizes the latent variables:
#'      diag( t(plscan_results$lx) %*% plscan_results$ly ) ## same as plsreg_results$d
#'      plscan_results$d
#'
#'      ## but this is an symmetric relationship:
#'      crossprod(plscan_results$lx) #orthogonal
#'      crossprod(plscan_results$ly) #orthogonal
#'      ## where the difference is in the maximization
#'      t(plscan_results$lx) %*% plscan_results$ly ## not orthogonal
#'  }
#'
#' @keywords multivariate, diagonalization, partial least squares


pls_can <- function(X, Y, center_X = TRUE, center_Y = TRUE, scale_X = TRUE, scale_Y = TRUE, components = 0, tol = .Machine$double.eps){


  X <- scale(X, center = center_X, scale = scale_X)
  if(center_X){
    X_center <- attributes(X)$`scaled:center`
  }else{
    X_center <- rep(0, ncol(X))
  }
  if(scale_X){
    X_scale <- attributes(X)$`scaled:scale`
  }else{
    X_scale <- rep(1, ncol(X))
  }

  Y <- scale(Y, center = center_Y, scale = scale_Y)
  if(center_Y){
    Y_center <- attributes(Y)$`scaled:center`
  }else{
    Y_center <- rep(0, ncol(Y))
  }
  if(scale_Y){
    Y_scale <- attributes(Y)$`scaled:scale`
  }else{
    Y_scale <- rep(1, ncol(Y))
  }


  gpls_can_results <- gpls_can(X = X, Y = Y, components = components, tol = tol)

  gpls_can_results$X_hat <- gpls_can_results$X_reconstructed * matrix(X_scale,nrow(X),ncol(X),byrow=T) + matrix(X_center,nrow(X),ncol(X),byrow=T)
  gpls_can_results$Y_hat <- gpls_can_results$Y_reconstructed * matrix(Y_scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y_center,nrow(Y),ncol(Y),byrow=T)

  gpls_can_results$X_hats <- array(NA,dim=c(nrow(X), ncol(X), length(gpls_can_results$d)))
  gpls_can_results$Y_hats <- array(NA,dim=c(nrow(Y), ncol(Y), length(gpls_can_results$d)))
  for(i in 1:length(gpls_can_results$d)){
    gpls_can_results$X_hats[,,i] <- gpls_can_results$X_reconstructeds[,,i] * matrix(X_scale,nrow(X),ncol(X),byrow=T) + matrix(X_center,nrow(X),ncol(X),byrow=T)
    gpls_can_results$Y_hats[,,i] <- gpls_can_results$Y_reconstructeds[,,i] * matrix(Y_scale,nrow(Y),ncol(Y),byrow=T) + matrix(Y_center,nrow(Y),ncol(Y),byrow=T)
  }

  return(gpls_can_results)

}
