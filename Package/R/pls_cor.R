#' @export
#'
#' @title Partial least squares correlation (PLSCOR)
#'
#' @description Computes partial least squares "correlation" between two data matrices by way of generalized PLS correlation
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
#' @return a list directly from \code{link{gpls_cor}}, which comes from \code{\link[GSVD]{gplssvd}}
#' \item{d.orig}{A vector containing the singular values of DAT above the tolerance threshold (based on eigenvalues).}
#' \item{l.orig}{A vector containing the eigen values of DAT above the tolerance threshold (\code{tol}).}
#' \item{tau}{A vector that contains the (original) explained variance per component (via eigenvalues: \code{$l.orig}.}
#' \item{d}{A vector of length \code{min(length(d.orig), k)} containing the retained singular values}
#' \item{l}{A vector of length \code{min(length(l.orig), k)} containing the retained eigen values}
#' \item{u}{Left (rows) singular vectors.}
#' \item{p}{Left (rows) generalized singular vectors.}
#' \item{fi}{Left (rows) component scores.}
#' \item{lx}{Latent variable scores for rows of \code{X}}
#' \item{v}{Right (columns) singular vectors.}
#' \item{q}{Right (columns) generalized singular vectors.}
#' \item{fj}{Right (columns) component scores.}
#' \item{ly}{Latent variable scores for rows of \code{Y}}
#'
#' @seealso \code{\link{cca}} \code{\link{rrr}} \code{\link{rda}} \code{\link[GSVD]{gplssvd}}
#'
#' @references
#' Abdi H., Eslami, A., Guillemot, V., & Beaton., D. (2018). Canonical correlation analysis (CCA). In R. Alhajj and J. Rokne (Eds.), \emph{Encyclopedia of Social Networks and Mining (2nd Edition).} New York: Springer Verlag.
#' Bookstein, F.L., 1994. Partial least squares: a dose–response model for measurement inthe behavioral and brain sciences. \emph{Psycoloquy} 5
#' Krishnan, A., Williams, L.J., McIntosh, A.R., & Abdi, H. (2011). Partial Least Squares (PLS) methods for neuroimaging: A tutorial and review. \emph{NeuroImage}, \bold{56}, 455-475.
#' McIntosh, A.R., Lobaugh, N.J., 2004. Partial least squares analysis of neuroimaging data: applications and advances. \emph{Neuroimage}, \bold{23}, S250–S263.
#' McIntosh, A.R., Bookstein, F., Haxby, J., Grady, C., 1996. Spatial pattern analysis offunctional brain images using partial least squares. \emph{Neuroimage}, \bold{3}, 143–157
#' Tucker, L., 1958. An inter-battery method of factor analysis. \emph{Psychometrika}, \bold{23}, 111–136
#'
#' @examples
#'
#'  data("wine", package = "GSVD")
#'  plscor_results <- pls_cor(wine$objective, wine$subjective)
#'
#'  \dontrun{
#'      ## like in RRR/RDA and CCA, PLSCOR maximizes the latent variables:
#'      t(plscor_results$lx) %*% plscor_results$ly ## a diagonal matrix
#'      diag( t(plscor_results$lx) %*% plscor_results$ly ) ## same as plscor_results$d
#'      plscor_results$d
#'  }
#'
#' @keywords multivariate, diagonalization, partial least squares

pls_cor <- function(X, Y, center_X = TRUE, center_Y = TRUE, scale_X = TRUE, scale_Y = TRUE, components = 0, tol = .Machine$double.eps){


  X <- scale(X, center = center_X, scale = scale_X)
  Y <- scale(Y, center = center_Y, scale = scale_Y)

  gpls_cor(X, Y,  components = components, tol = tol)

}
