#' @export
#'
#' @title Partial least squares-correspondence analysis "regression decomposition" (PLSCAREG)
#'
#' @description Computes partial least squares-correspondence analysis "regression decomposition" between two data matrices by way of generalized PLS regression decomposition
#'  Each data matrix is assumed to be suitable for correspondence analysis and preprocessed (e.g., categorical through disjunctive coding; ordinal through thermometer coding; continuous through Escofier coding)
#'
#' @param X Data matrix with \emph{I} rows and \emph{J} columns
#' @param Y Data matrix with \emph{I} rows and \emph{K} columns
#' @param components The number of components to return. If < 1 then the maximum components will be returned. Default = 0.
#' @param tol default is .Machine$double.eps. A parameter to pass through to \code{\link[GSVD]{gplssvd}}; eliminates effectively zero, negative, or imaginary singular values (thus components).
#'
#' @return A list of outputs based on \code{\link{gpls_reg}}
#' \item{d}{A vector containing the singular values from each iteration.}
#' \item{u}{Left (rows) singular vectors.}
#' \item{v}{Right (columns) singular vectors. In PLSREG sometimes called "weight matrix".}
#' \item{lx}{Latent variable scores for rows of \code{X}}
#' \item{ly}{Latent variable scores for rows of \code{Y}}
#' \item{p}{Left (rows) generalized singular vectors.}
#' \item{q}{Right (columns) generalized singular vectors.}
#' \item{fi}{Left (rows) component scores.}
#' \item{fj}{Right (columns) component scores.}
#' \item{tx}{"Latent vectors": A normed version of \code{lx} for use in rebuilding \code{X} data}
#' \item{u_hat}{"Loading matrix": A "predicted" version of \code{u} for use in rebuilding \code{X} data}
#' \item{betas}{"Regression weights": Akin to betas for use in rebuilding \code{Y}}
#' \item{X_reconstructeds}{A version of \code{X} reconstructed for each iteration (i.e., latent variable/component)}
#' \item{Y_reconstructeds}{A version of \code{Y} reconstructed for each iteration (i.e., latent variable/component)}
#' \item{X_residuals}{The residualized (i.e., \code{X - X_reconstructeds}) version of \code{X} for each iteration (i.e., latent variable/component)}
#' \item{Y_residuals}{The residualized (i.e., \code{Y - Y_reconstructeds}) version of \code{Y} for each iteration (i.e., latent variable/component)}
#' \item{r2_x}{Proporition of explained variance from \code{X} to each latent variable/component.}
#' \item{r2_y}{Proporition of explained variance from \code{Y} to each latent variable/component.}
#' \item{Y_reconstructed}{A version of \code{Y} reconstructed from all iterations (i.e., latent variables/components); see \code{components}.}
#' \item{Y_residual}{The residualized (i.e., \code{Y - Y_reconstructed} from all iterations (i.e., latent variables/components); see \code{components}.}
#' \item{Y_hat}{The version of \code{Y_reconstructed} in the original coordinates (expected & marginal values) as \code{Y}.}
#' \item{X_hats}{The versions of \code{X_reconstructeds} in the original coordinates (expected & marginal values) as \code{X} per iteration (i.e., latent variable/component).}
#' \item{Y_hats}{The versions of \code{Y_reconstructeds} in the original coordinates (expected & marginal values) as \code{Y} per iteration (i.e., latent variable/component).}
#'
#' @seealso \code{\link{ca_preproc}} \code{\link{make_data_disjunctive}} \code{\link{thermometer_coding}} \code{\link{escofier_coding}} \code{\link{pls_reg}} \code{\link{plsca_cor}} \code{\link{plsca_can}} \code{\link{gpls_reg}} \code{\link[GSVD]{gplssvd}}
#'
#' @references
#' Beaton, D., Dunlop, J., & Abdi, H. (2016). Partial least squares correspondence analysis: A framework to simultaneously analyze behavioral and genetic data. \emph{Psychological methods}, \bold{21} (4), 621.
#' Beaton, D., ADNI, Saporta, G., Abdi, H. (2019). A generalization of partial least squares regression and correspondence analysis for categorical and mixed data: An application with the ADNI data. \emph{bioRxiv}, 598888.
#'
#' @examples
#'
#'  data("snps.druguse", package = "GSVD")
#'  X <- make_data_disjunctive(snps.druguse$DATA1)
#'  Y <- make_data_disjunctive(snps.druguse$DATA2)
#'
#'  plscareg_results <- plsca_reg(X, Y)
#'
#'  \dontrun{
#'      ## like in PLSREG, PLSCAREG maximizes the latent variables:
#'      diag( t(plscareg_results$lx) %*% plscareg_results$ly ) ## same as plscareg_results$d
#'      plscareg_results$d
#'
#'      ## but this is an asymmetric relationship:
#'      crossprod(plscareg_results$lx) #orthogonal
#'      t(plscareg_results$lx) %*% plscareg_results$ly ## a triangular matrix
#'  }
#'
#' @keywords multivariate, diagonalization, partial least squares

plsca_reg <- function(X, Y, components = 0, tol = .Machine$double.eps){

  X_ca_preproc <- ca_preproc(X)
  Y_ca_preproc <- ca_preproc(Y)


  gpls_reg_results <- gpls_reg(X_ca_preproc$Z, Y_ca_preproc$Z,
                               XLW = diag(1/X_ca_preproc$m), XRW = diag(1/X_ca_preproc$w),
                               YLW = diag(1/Y_ca_preproc$m), YRW = diag(1/Y_ca_preproc$w),
                               components = components, tol = tol)

    ### not exactly the same residual as in PLS REG; should harmonize or check this.
  gpls_reg_results$Y_residual <- ((Y_ca_preproc$Z - gpls_reg_results$Y_reconstructed) + Y_ca_preproc$E) * sum(Y)
  gpls_reg_results$Y_hat <- (gpls_reg_results$Y_reconstructed + Y_ca_preproc$E) * sum(Y)

  gpls_reg_results$X_hats <- array(NA,dim=c(nrow(X), ncol(X), length(gpls_reg_results$d)))
  gpls_reg_results$Y_hats <- array(NA,dim=c(nrow(Y), ncol(Y), length(gpls_reg_results$d)))

  for(i in 1:length(gpls_reg_results$d)){
     gpls_reg_results$X_hats[,,i] <- (gpls_reg_results$X_reconstructeds[,,i] + X_ca_preproc$E) * sum(X)
     gpls_reg_results$Y_hats[,,i] <- (gpls_reg_results$Y_reconstructeds[,,i] + Y_ca_preproc$E) * sum(Y)
  }

  rownames(gpls_reg_results$X_hats) <- rownames(X)
  colnames(gpls_reg_results$X_hats) <- colnames(X)

  rownames(gpls_reg_results$Y_hat) <- rownames(gpls_reg_results$Y_residual) <- rownames(gpls_reg_results$Y_hats) <- rownames(Y)
  colnames(gpls_reg_results$Y_hat) <- colnames(gpls_reg_results$Y_residual) <- colnames(gpls_reg_results$Y_hats) <- colnames(Y)

  return(gpls_reg_results)


}
