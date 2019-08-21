#' @export
#'
#' @title Partial least squares-correspondence analysis correlation (PLSCACOR)
#'
#' @description Computes partial least squares-correspondence analysis "correlation" between two data matrices by way of generalized PLS correlation
#'  Each data matrix is assumed to be suitable for correspondence analysis and preprocessed (e.g., categorical through disjunctive coding; ordinal through thermometer coding; continuous through Escofier coding)
#'
#' @param X Data matrix with \emph{I} rows and \emph{J} columns
#' @param Y Data matrix with \emph{I} rows and \emph{K} columns
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
#' @seealso \code{\link{ca_preproc}} \code{\link{make_data_disjunctive}} \code{\link{thermometer_coding}} \code{\link{escofier_coding}} \code{\link{pls_cor}} \code{\link{plsca_reg}} \code{\link{plsca_can}} \code{\link[GSVD]{gplssvd}}
#'
#' @references
#' Beaton, D., Dunlop, J., & Abdi, H. (2016). Partial least squares correspondence analysis: A framework to simultaneously analyze behavioral and genetic data. \emph{Psychological methods}, \bold{21} (4), 621.
#'
#' @examples
#'
#'  data("snps.druguse", package = "GSVD")
#'  X <- make_data_disjunctive(snps.druguse$DATA1)
#'  Y <- make_data_disjunctive(snps.druguse$DATA2)
#'
#'  plscacor_results <- plsca_cor(X, Y)
#'
#'  \dontrun{
#'      ## like in PLSCOR, PLSCACOR maximizes the latent variables:
#'      t(plscacor_results$lx) %*% plscacor_results$ly ## a diagonal matrix
#'      diag( t(plscacor_results$lx) %*% plscacor_results$ly ) ## same as plscacor_results$d
#'      plscacor_results$d
#'  }
#'
#' @keywords multivariate, diagonalization, partial least squares

plsca_cor <- function(X, Y, components = 0, tol = .Machine$double.eps){

  X_ca_preproc <- ca_preproc(X)
  Y_ca_preproc <- ca_preproc(Y)

  gpls_cor(X_ca_preproc$Z, Y_ca_preproc$Z,
          XLW = diag(1/X_ca_preproc$m), XRW = diag(1/X_ca_preproc$w),
          YLW = diag(1/Y_ca_preproc$m), YRW = diag(1/Y_ca_preproc$w),
          components = components, tol = tol)


}
