#' @export
#'
#' @title Reduced rank regression (see also redundancy analysis)
#'
#' @description Computes reduced rank regression between two data matrices by way of generalized PLS correlation
#'
#' @param X Data matrix with \emph{I} rows and \emph{J} columns
#' @param Y Data matrix with \emph{I} rows and \emph{K} columns
#' @param center_X For the \code{X} matrix: A parameter to pass through to \code{center} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{X}.
#' @param center_Y For the \code{Y} matrix: A parameter to pass through to \code{center} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{Y}.
#' @param scale_X For the \code{X} matrix: A parameter to pass through to \code{scale} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{X}.
#' @param scale_Y For the \code{Y} matrix: A parameter to pass through to \code{scale} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{Y}.
#' @param components The number of components to return. If < 1 then the maximum components will be returned. Default = 0.
#' @param tol default is .Machine$double.eps. A parameter to pass through to \code{\link[GSVD]{gplssvd}}; eliminates singular values that are effectively zero (and thus drops null components).
#'
#'
#' @return a list directly from \code{link{gpls_cor}}, which comes from \code{\link[GSVD]{gplssvd}}, with one additional item
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
#' \item{beta_matrix}{A matrix of beta estimates}
#'
#' @seealso \code{\link{pls_cor}} \code{\link{rda}} \code{\link{cca}} \code{\link[GSVD]{gplssvd}}
#'
#' @note
#' "Reduced rank regression", also known as "redundancy analysis" (amongst other names) optimizes for latent variables in a way somewhere between PLS correlation and CCA.
#' See the references for further discussions. See also an older implementation in S: \url{http://ftp.uni-bayreuth.de/math/statlib/S/rrr.s}
#'
#' Additionally, this implementation makes use of generalized/psuedo inverse for computations, which allows RRR/RDA to run even when data
#' do not meet assumptions typically expected for standard RRR/RDA, e.g., collinear, rank deficient.
#'
#' @references
#' Abdi H., Eslami, A., Guillemot, V., & Beaton., D. (2018). Canonical correlation analysis (CCA). In R. Alhajj and J. Rokne (Eds.), \emph{Encyclopedia of Social Networks and Mining (2nd Edition).} New York: Springer Verlag.
#'
#'
#' @examples
#'
#'  data("wine", package = "GSVD")
#'  rrr_results <- rrr(wine$objective, wine$subjective)
#'
#'  \dontrun{
#'      ## like in CCA and PLSCOR, RRR maximizes the latent variables:
#'      t(rrr_results$lx) %*% rrr_results$ly ## a diagonal matrix
#'      diag( t(rrr_results$lx) %*% rrr_results$ly ) ## same as rrr_results$d
#'      rrr_results$d
#'  }
#'
#' @keywords multivariate, diagonalization, reduced rank regression, redundancy analysis
#' @importFrom MASS ginv


rrr <- function(X, Y, center_X = TRUE, center_Y = TRUE, scale_X = TRUE, scale_Y = TRUE, components = 0, tol = .Machine$double.eps){

  X <- scale(X, center = center_X, scale = scale_X)
  Y <- scale(Y, center = center_Y, scale = scale_Y)

  ## a bit of a trick
  rrr_res <- gpls_cor(
    X = t(MASS::ginv(X)), Y = Y,
    XRW=crossprod(X),
    components = components, tol = tol
  )

  ## I should send back "Betas" here per http://ftp.uni-bayreuth.de/math/statlib/S/rrr.s
    ### and  I should consider including those into CCA and PLSC
    ### and also for the other side, and also consider the rebuilding procedure here as is done in gpls_can and gpls_reg
  rrr_res$beta_matrix <- (rrr_res$p %*% diag(rrr_res$d))

  rrr_res
}


#' @export
#'
#' @title Redundancy analysis (see also reduced rank regression)
#'
#' @description Computes redundancy analysis between two data matrices by way of generalized PLS correlation
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
#' @return a list directly from \code{link{gpls_cor}}, which comes from \code{\link[GSVD]{gplssvd}}, with one additional item
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
#' \item{beta_matrix}{A matrix of beta estimates}
#'
#' @seealso \code{\link{pls_cor}} \code{\link{rda}} \code{\link{cca}} \code{\link[GSVD]{gplssvd}}
#'
#' @references
#' Abdi H., Eslami, A., Guillemot, V., & Beaton., D. (2018). Canonical correlation analysis (CCA). In R. Alhajj and J. Rokne (Eds.), \emph{Encyclopedia of Social Networks and Mining (2nd Edition).} New York: Springer Verlag.
#' Rao C. R. (1964) The use and interpretation of principalcomponent analysis in applied research. \emph{Sankhyā}, A26:329–358
#' Van Den Wollenberg A. L. (1977) Redundancy analysis analternative  for  canonical  correlation  analysis. \emph{Psychometrika}, 42:207–219
#'
#' @note
#' "Redundancy analysis", also known as "reduced rank regression" (amongst other names) optimizes for latent variables in a way somewhere between PLS correlation and CCA.
#' See the references for further discussions. See also an older implementation in S: \url{http://ftp.uni-bayreuth.de/math/statlib/S/rrr.s}
#'
#' @examples
#'
#'  data("wine", package = "GSVD")
#'  rda_results <- rda(wine$objective, wine$subjective)
#'
#'  \dontrun{
#'      ## like in CCA and PLSCOR, RDA maximizes the latent variables:
#'      t(rda_results$lx) %*% rda_results$ly ## a diagonal matrix
#'      diag( t(rda_results$lx) %*% rda_results$ly ) ## same as rda_results$d
#'      rda_results$d
#'  }
#'
#'
#' @keywords multivariate, diagonalization, reduced rank regression, redundancy analysis

rda <- function(X, Y, center_X = TRUE, center_Y = TRUE, scale_X = TRUE, scale_Y = TRUE, components = 0, tol = .Machine$double.eps){

  rrr(X = X, Y = Y, center_X = center_X, center_Y = center_Y, scale_X = scale_X, scale_Y = scale_Y, components = components, tol = tol)

}
