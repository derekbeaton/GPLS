#' @export
#'
#' @title Generalized partial least squares "canonical decomposition" (GPLSCAN)
#'
#' @description Computes generalized partial least squares "canonical decomposition" between two data matrices.
#' GPLSCAN allows for the use of left (row) and right (column) weights for each data matrix.
#'
#' @param X Data matrix with \emph{I} rows and \emph{J} columns
#' @param Y Data matrix with \emph{I} rows and \emph{K} columns
#' @param XLW An \emph{I} by \emph{I} matrix of row weights for \code{X}. Default is \code{diag(nrow(X))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param YLW An \emph{I} by \emph{I} matrix of row weights for \code{Y}. Default is \code{diag(nrow(Y))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param XRW A \emph{J} by \emph{J} matrix of row weights for \code{X}. Default is \code{diag(ncol(X))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param YRW A \emph{K} by \emph{K} matrix of row weights for \code{Y}. Default is \code{diag(ncol(Y))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param components The number of components to return. If < 1 then the maximum components will be returned. Default = 0.
#' @param tol default is .Machine$double.eps. A parameter to pass through to \code{\link[GSVD]{gplssvd}}; eliminates effectively zero, negative, or imaginary singular values (thus components).
#'
#' @return A list of outputs
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
#'
#' @seealso \code{\link{gpls_reg}} \code{\link{gpls_cor}} \code{\link{pls_can}} \code{\link{plsca_can}} \code{\link[GSVD]{gplssvd}}
#'
#' @references
#' Beaton, D., ADNI, Saporta, G., Abdi, H. (2019). A generalization of partial least squares regression and correspondence analysis for categorical and mixed data: An application with the ADNI data. \emph{bioRxiv}, 598888.
#'
#' @examples
#'
#'  \dontrun{
#'  library(GSVD)
#'  data("wine", package = "GSVD")
#'  X <- scale(wine$objective)
#'  Y <- scale(wine$subjective)
#'
#'
#'  ##partial least squares "regression decomposition"
#'  #### the first latent variable from reg & cor & can are identical in all PLSs.
#'  gplscan_pls_optimization <- gpls_can( X, Y)
#'
#'  ## "partial least squares regression"
#'  ### but with the optimization per latent variable of CCA
#'  #### because of optimization, this ends up identical to
#'  #### cca(X, Y, center_X = F, center_Y = F, scale_X = F, scale_Y = F)
#'  gplscan_cca_optimization <- gpls_can( t(MASS::ginv(X)), t(MASS::ginv(Y)),
#'       XRW = crossprod(X), YRW = crossprod(Y))
#'
#'  ## "partial least squares regression"
#'  ### but with the optimization per latent variable of RRR/RDA
#'  #### because of optimization, this ends up NEARLY identical to
#'  #### rrr(X, Y, center_X = F, center_Y = F, scale_X = F, scale_Y = F)
#'  #### or rda(X, Y, center_X = F, center_Y = F, scale_X = F, scale_Y = F)
#'  ##### the only difference between common items is $ly
#'  gplscan_rrr_optimization <- gpls_can( t(MASS::ginv(X)), Y, XRW = crossprod(X))
#'
#'  rm(X)
#'  rm(Y)
#'
#'  ## standard partial least squares-correspondence analysis "canonical decomposition"
#'  #### the first latent variable from reg & cor & can are identical in all PLSs.
#'  data("snps.druguse", package = "GSVD")
#'  X <- make_data_disjunctive(snps.druguse$DATA1)
#'  Y <- make_data_disjunctive(snps.druguse$DATA2)
#'
#'  X_ca_preproc <- ca_preproc(X)
#'  Y_ca_preproc <- ca_preproc(Y)
#'
#'  gplscan_plsca <- gpls_can( X = X_ca_preproc$Z, Y = Y_ca_preproc$Z,
#'      XLW = diag(1/X_ca_preproc$m), YLW = diag(1/Y_ca_preproc$m),
#'      XRW = diag(1/X_ca_preproc$w), YRW = diag(1/Y_ca_preproc$w)
#'  )
#'  }
#'
#' @keywords multivariate, diagonalization, partial least squares
#' @importFrom MASS ginv
#' @importFrom GSVD invsqrt_psd_matrix sqrt_psd_matrix

gpls_can <- function(X, Y,
                     XLW = diag(nrow(X)), YLW = diag(nrow(Y)),
                     XRW = diag(ncol(X)), YRW = diag(ncol(Y)),
                     components = 0, tol = .Machine$double.eps){


  X_gsvd <- gsvd(X, XLW, XRW)
    X_rank <- length(X_gsvd$d)
    X_trace <- sum(X_gsvd$d^2)
  rm(X_gsvd)

  Y_gsvd <- gsvd(Y, YLW, YRW)
    Y_rank <- length(Y_gsvd$d)
    Y_trace <- sum(Y_gsvd$d^2)
  rm(Y_gsvd)

  stopped_early <- F

  if((components > X_rank) | (components > Y_rank) | (components < 1)){
    components <- min(X_rank, Y_rank)
  }

  tx <- lx <- matrix(NA,nrow(X),components)
  ty <- ly <- matrix(NA,nrow(Y),components)

  u_hat <- fi <- p <- u <- matrix(NA,ncol(X),components)
  v_hat <- fj <- q <- v <- matrix(NA,ncol(Y),components)

  r2_x_cumulative <- r2_y_cumulative <- d <- rep(NA, components)

  X_reconstructeds <- X_residuals <- array(NA,dim=c(nrow(X),ncol(X),components))
  Y_reconstructeds <- Y_residuals <- array(NA,dim=c(nrow(Y),ncol(Y),components))

  X_deflate <- X
  Y_deflate <- Y

  for(i in 1:components){

    gplssvd_results <- gpls_cor(X_deflate, Y_deflate,
                                XRW = XRW, YRW = YRW,
                                XLW = XLW, YLW = YLW,
                                components = 1, tol = tol)

    u[,i] <- gplssvd_results$u
    p[,i] <- gplssvd_results$p
    fi[,i] <- gplssvd_results$fi
    v[,i] <- gplssvd_results$v
    q[,i] <- gplssvd_results$q
    fj[,i] <- gplssvd_results$fj
    d[i] <- gplssvd_results$d
    lx[,i] <- gplssvd_results$lx
    ly[,i] <- gplssvd_results$ly

    tx[,i] <- lx[,i] / sqrt(sum(lx[,i]^2))
    u_hat[,i] <- t(tx[,i]) %*% ( GSVD::sqrt_psd_matrix(XLW) %*% X_deflate %*% GSVD::sqrt_psd_matrix(XRW) )

    ty[,i] <- ly[,i] / sqrt(sum(ly[,i]^2))
    v_hat[,i] <- t(ty[,i]) %*% ( GSVD::sqrt_psd_matrix(YLW) %*% Y_deflate %*% GSVD::sqrt_psd_matrix(YRW) )

    X_reconstructeds[,,i] <- GSVD::invsqrt_psd_matrix(XLW) %*% (tx[,i] %o% u_hat[,i]) %*% GSVD::invsqrt_psd_matrix(XRW)
      X_reconstructeds[abs(X_reconstructeds) < tol] <- 0
    Y_reconstructeds[,,i] <- GSVD::invsqrt_psd_matrix(YLW) %*% (ty[,i] %o% v_hat[,i]) %*% GSVD::invsqrt_psd_matrix(YRW)
      Y_reconstructeds[abs(Y_reconstructeds) < tol] <- 0

    X_residuals[,,i] <- (X_deflate - X_reconstructeds[,,i])
      X_residuals[abs(X_residuals) < tol] <- 0
    Y_residuals[,,i] <- (Y_deflate - Y_reconstructeds[,,i])
      Y_residuals[abs(Y_residuals) < tol] <- 0

    X_deflate <- X_residuals[,,i]
    Y_deflate <- Y_residuals[,,i]


    r2_x_cumulative[i] <- (X_trace-sum( ( GSVD::sqrt_psd_matrix(XLW) %*%  X_deflate %*% GSVD::sqrt_psd_matrix(XRW) ) ^2)) / X_trace
    r2_y_cumulative[i] <- (Y_trace-sum( ( GSVD::sqrt_psd_matrix(YLW) %*%  Y_deflate %*% GSVD::sqrt_psd_matrix(YRW) ) ^2)) / Y_trace


    if( (sum(Y_deflate^2) < tol) & (i < components) ){

      stopped_early <- T
      warning("gpls_can: Y is fully deflated. Stopping early.")

    }
    if( (sum(X_deflate^2) < tol) & (i < components) ){

      stopped_early <- T
      warning("gpls_can: X is fully deflated. Stopping early.")

    }

    if(stopped_early){
      break
    }
  }

  if(stopped_early){
    u <- u[,1:i]
    p <- p[,1:i]
    fi <- fi[,1:i]
    v <- v[,1:i]
    q <- q[,1:i]
    fj <- fj[,1:i]
    d <- d[1:i]
    lx <- lx[,1:i]
    ly <- ly[,1:i]

    u_hat <- u_hat[,1:i]
    v_hat <- v_hat[,1:i]

    tx = tx[,1:i]
    ty = ty[,1:i]

    X_reconstructeds <- X_reconstructeds[,,1:i]
    Y_reconstructeds <- Y_reconstructeds[,,1:i]
    X_residuals <- X_residuals[,,1:i]
    Y_residuals <- Y_residuals[,,1:i]

    r2_x_cumulative <- r2_x_cumulative[1:i]
    r2_y_cumulative <- r2_y_cumulative[1:i]
  }


  X_reconstructed <- GSVD::invsqrt_psd_matrix(XLW) %*% (tx %*% t(u_hat)) %*% GSVD::invsqrt_psd_matrix(XRW)
    X_reconstructed[abs(X_reconstructed) < tol] <- 0
  X_residual <- X - X_reconstructed

  Y_reconstructed <- GSVD::invsqrt_psd_matrix(YLW) %*% (ty %*% t(v_hat)) %*% GSVD::invsqrt_psd_matrix(YRW)
    Y_reconstructed[abs(Y_reconstructed) < tol] <- 0
  Y_residual <- Y - Y_reconstructed

  rownames(tx) <- rownames(lx) <- rownames(X_reconstructed) <- rownames(X_residual) <- rownames(X_reconstructeds) <- rownames(X_residuals) <- rownames(X)
  rownames(fi) <- rownames(u) <- rownames(p) <- rownames(u_hat) <- colnames(X_reconstructed) <- colnames(X_residual) <- colnames(X_reconstructeds) <- colnames(X_residuals) <- colnames(X)
  rownames(ty) <- rownames(ly) <- rownames(Y_reconstructed) <- rownames(Y_residual) <- rownames(Y_reconstructeds) <- rownames(Y_residuals) <- rownames(Y)
  rownames(fj) <- rownames(v) <- rownames(q) <- rownames(v_hat) <- colnames(Y_reconstructed) <- colnames(Y_residual) <- colnames(Y_reconstructeds) <- colnames(Y_residuals) <- colnames(Y)

  return( list(
    d = d, u = u, v = v, lx = lx, ly = ly,
    p = p, q = q, fi = fi, fj = fj,
    tx = tx, ty = ty,
    u_hat = u_hat, v_hat = v_hat,
    X_reconstructeds = X_reconstructeds, X_residuals = X_residuals,
    Y_reconstructeds = Y_reconstructeds, Y_residuals = Y_residuals,
    r2_x = diff(c(0,r2_x_cumulative)), r2_y = diff(c(0,r2_y_cumulative)),
    X_reconstructed = X_reconstructed, X_residual = X_residual,
    Y_reconstructed = Y_reconstructed, Y_residual = Y_residual
  ) )


}
