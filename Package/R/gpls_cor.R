#' @export
#'
#' @title Generalized partial least squares "correlation decomposition" (GPLSCOR)
#'
#' @description Computes generalized partial least squares "correlation decomposition" between two data matrices.
#' GPLSCOR allows for the use of left (row) and right (column) weights for each data matrix.
#'
#' @param X Data matrix with \emph{I} rows and \emph{J} columns
#' @param Y Data matrix with \emph{I} rows and \emph{K} columns
#' @param XLW An \emph{I} by \emph{I} matrix of row weights for \code{X}. Default is \code{diag(nrow(X))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param YLW An \emph{I} by \emph{I} matrix of row weights for \code{Y}. Default is \code{diag(nrow(Y))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param XRW A \emph{J} by \emph{J} matrix of row weights for \code{X}. Default is \code{diag(ncol(X))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param YRW A \emph{K} by \emph{K} matrix of row weights for \code{Y}. Default is \code{diag(ncol(Y))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param components The number of components to return. If < 1 then the maximum components will be returned. Default = 0.
#' @param tol default is .Machine$double.eps. A parameter to pass through to \code{\link[GSVD]{gplssvd}}; eliminates singular values that are effectively zero (and thus drops null components).
#'
#'
#' @return a list of outputs; see also \code{\link[GSVD]{gplssvd}}
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
#' @seealso \code{\link{gpls_reg}} \code{\link{gpls_can}} \code{\link{pls_cor}} \code{\link{plsca_cor}} \code{\link{cca}} \code{\link{rrr}} \code{\link{rda}} \code{\link[GSVD]{gplssvd}}
#'
#' @references
#' Beaton, D., ADNI, Saporta, G., Abdi, H. (2019). A generalization of partial least squares regression and correspondence analysis for categorical and mixed data: An application with the ADNI data. \emph{bioRxiv}, 598888.
#' Beaton, D., Dunlop, J., & Abdi, H. (2016). Partial least squares correspondence analysis: A framework to simultaneously analyze behavioral and genetic data. \emph{Psychological methods}, \bold{21} (4), 621.
#'
#' @examples
#'
#'  \dontrun{
#'  library(GSVD)
#'  data("wine", package = "GSVD")
#'  X <- scale(wine$objective)
#'  Y <- scale(wine$subjective)
#'
#'  ## partial least squares "correlation decomposition"
#'  #### the first latent variable from reg & cor & can are identical in all PLSs.
#'  #### this is pls_cor(X, Y, center_X = F, center_Y = F, scale_X = F, scale_Y = F)
#'  gplscor_pls_optimization <- gpls_cor(X, Y)
#'
#'  ## partial least squares "correlation decomposition"
#'  ### but with the optimization per latent variable of CCA
#'  #### this is cca(X, Y, center_X = F, center_Y = F, scale_X = F, scale_Y = F)
#'  gplscor_cca_optimization <- gpls_cor( t(MASS::ginv(X)), t(MASS::ginv(Y)),
#'       XRW = crossprod(X), YRW = crossprod(Y))
#'
#'  ## partial least squares "correlation decomposition"
#'  ### but with the optimization per latent variable of RRR/RDA
#'  #### this is rrr(X, Y, center_X = F, center_Y = F, scale_X = F, scale_Y = F) or
#'  #### rda(X, Y, center_X = F, center_Y = F, scale_X = F, scale_Y = F)
#'  gplscor_rrr_optimization <- gpls_cor( t(MASS::ginv(X)), Y, XRW = crossprod(X))
#'
#'
#'  rm(X)
#'  rm(Y)
#'
#'  ## partial least squares-correspondence analysis "correlation decomposition"
#'  #### the first latent variable from reg & cor & can are identical in all PLSs.
#'  data("snps.druguse", package = "GSVD")
#'  X <- make_data_disjunctive(snps.druguse$DATA1)
#'  Y <- make_data_disjunctive(snps.druguse$DATA2)
#'
#'  X_ca_preproc <- ca_preproc(X)
#'  Y_ca_preproc <- ca_preproc(Y)
#'
#'  #### this is plsca_cor(X, Y)
#'  gplscor_plsca <- gpls_cor( X = X_ca_preproc$Z, Y = Y_ca_preproc$Z,
#'      XLW = diag(1/X_ca_preproc$m), YLW = diag(1/Y_ca_preproc$m),
#'      XRW = diag(1/X_ca_preproc$w), YRW = diag(1/Y_ca_preproc$w)
#'  )
#'
#'  }
#'
#' @keywords multivariate, diagonalization, partial least squares
#' @importFrom MASS ginv

gpls_cor <- function(X, Y,
                  XLW = diag(nrow(X)), YLW = diag(nrow(Y)),
                  XRW = diag(ncol(X)), YRW = diag(ncol(Y)),
                  components = 0, tol = .Machine$double.eps){

    ## that's it.
    gplssvd(X, Y, XLW = XLW, YLW = YLW, XRW = XRW, YRW = YRW, k = components, tol = tol)

}
