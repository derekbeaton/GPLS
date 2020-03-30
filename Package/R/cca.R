#' @export
#'
#' @title Canonical correlation analysis (CCA)
#'
#' @description Computes canonical correlation analysis between two data matrices by way of generalized PLS correlation
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
#' @seealso \code{\link{pls_cor}} \code{\link{rrr}} \code{\link{rda}} \code{\link{gpls_cor}} \code{\link[GSVD]{gplssvd}}
#'
#' @note
#' Additionally, this implementation makes use of generalized/psuedo inverse for computations, which allows CCA to run even when data
#' do not meet assumptions typically expected for standard CCA, e.g., collinear, rank deficient.
#'
#' @references
#' Abdi H., Eslami, A., Guillemot, V., & Beaton., D. (2018). Canonical correlation analysis (CCA). In R. Alhajj and J. Rokne (Eds.), \emph{Encyclopedia of Social Networks and Mining (2nd Edition).} New York: Springer Verlag.
#'
#'
#' @examples
#'
#'  data("wine", package = "GSVD")
#'  cca_results <- cca(wine$objective, wine$subjective)
#'
#'  \dontrun{
#'      ## like in PLSCOR and RRR/RDA, CCA maximizes the latent variables:
#'      t(cca_results$lx) %*% cca_results$ly ## a diagonal matrix
#'      diag( t(cca_results$lx) %*% cca_results$ly ) ## same as cca_results$d
#'      cca_results$d
#'
#'      ## where the singular values are the "canonical correlations"
#'          #### see \code{$cor} from \code{cancor}.
#'      cca_results$d
#'
#'      ## additionally, generalized singular vectors are the "coefficients"
#'          #### see \code{$xcoef} and \code{$ycoef} from \code{cancor}.
#'      cca_results$p
#'      cca_results$q
#'
#'  }
#'
#'
#' @keywords multivariate, diagonalization, canonical correlation
#' @importFrom MASS ginv

cca <- function(X, Y, center_X = TRUE, center_Y = TRUE, scale_X = TRUE, scale_Y = TRUE, components = 0, tol = .Machine$double.eps){


  X <- scale(X, center = center_X, scale = scale_X)
  Y <- scale(Y, center = center_Y, scale = scale_Y)

  ## a bit of a trick
  gpls_cor(
    X = MASS::ginv(X), Y = MASS::ginv(Y),
    XRW=crossprod(X), YRW=crossprod(Y),
    components = components, tol = tol
  )


}
