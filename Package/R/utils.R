#' GPLS
#'
#' @description A small package dedicated to multiple generalizations of partial least squares.
#' The first set of generalizations are optimizations such as canonical correlation (correlation) reduced rank regression a.k.a. redundancy analysis (regression), and partial least squares (covariance).
#' The second set of generalizations are that these optimizations can be performed through multiple PLS-based algorithms such as "PLSSVD" which is a single pass of the generalized singular value decomposition (GSVD), or "canonical partial least squares" which uses a symmetric deflation---of both X and Y---with an iterative use of the GSVD, and "partial least squares regression" which uses an asymmetric deflation---of only X with a prediction of Y from X---with an iterative use of the GSVD.
#' The third set of generalizations are that optimizations and algorithms can be used under assumptions of Correspondence Analysis (CA) which allows us to use PLS algorithms for continuous, categorical, ordinal, or mixed data types, all generally under the assumptions of normal (PCA) or Chi-squared (CA).
#' Effectively, GPLS allows for virtually any data type, any optmization (weights), and various algorithms all unified under the partial least squares and correspondence analysis nomenclatures and approaches.
#' @seealso \code{\link[GSVD]{gplssvd}}, \code{\link[GSVD]{gsvd}}, \code{\link{cca}}, \code{\link{rrr}}, \code{\link{pls_cor}}, \code{\link{pls_can}}, \code{\link{pls_reg}}, \code{\link{plsca_cor}}, \code{\link{plsca_can}}, \code{\link{plsca_reg}}
#'
#' @references
#' The detailed explanation of the generalizations can be found in:
#' #' Beaton, D., ADNI, Saporta, G., Abdi, H. (2019). A generalization of partial least squares regression and correspondence analysis for categorical and mixed data: An application with the ADNI data. \emph{bioRxiv}, 598888.
#'
#' Please also see the long history of partial least squares:
#' Abdi H., Eslami, A., Guillemot, V., & Beaton., D. (2018). Canonical correlation analysis (CCA). In R. Alhajj and J. Rokne (Eds.), \emph{Encyclopedia of Social Networks and Mining (2nd Edition).} New York: Springer Verlag.
#' Beaton, D., Dunlop, J., & Abdi, H. (2016). Partial least squares correspondence analysis: A framework to simultaneously analyze behavioral and genetic data. \emph{Psychological methods}, \bold{21} (4), 621.
#' Bookstein, F.L., 1994. Partial least squares: a dose–response model for measurement inthe behavioral and brain sciences. \emph{Psycoloquy} 5
#' Krishnan, A., Williams, L.J., McIntosh, A.R., & Abdi, H. (2011). Partial Least Squares (PLS) methods for neuroimaging: A tutorial and review. \emph{NeuroImage}, \bold{56}, 455-475.
#' McIntosh, A.R., Lobaugh, N.J., 2004. Partial least squares analysis of neuroimaging data: applications and advances. \emph{Neuroimage}, \bold{23}, S250–S263.
#' McIntosh, A.R., Bookstein, F., Haxby, J., Grady, C., 1996. Spatial pattern analysis offunctional brain images using partial least squares. \emph{Neuroimage}, \bold{3}, 143–157
#' Tenenhaus, M. (1998). La Regression PLS. Theorie et Pratique. \emph{Editions TECHNIP}, Paris.
#' Tucker, L., 1958. An inter-battery method of factor analysis. \emph{Psychometrika}, \bold{23}, 111–136
#'
#'
#' @docType package
#' @name GPLS
#' @aliases gpls gpls-package GPLS-package
#'
#' @import GSVD
#'
"_PACKAGE"


#' @export
#'
#' @title Correspondence Analysis (CA) preprocessing
#'
#' @description Standard approach to data preprocessing in correspondence analysis, i.e., deviations from independence just as in Chi-squared.
#'     Returns multiple matrices and vectors required by correspondence analysis
#'
#' @param DATA a matrix for preprocessing
#'
#' @return A list with five elements:
#' \item{O}{Matrix: the "observed" matrix of probabilities}
#' \item{m}{Vector: the row marginal probabilities of the "observed" matrix}
#' \item{w}{Vector: the column marginal probabilities of the "observed" matrix}
#' \item{E}{Matrix: the "expected" probabilities matrix}
#' \item{Z}{Matirx: the "deviations" matrix}

ca_preproc <- function(DATA){

  O <- DATA/sum(DATA)
  m <- rowSums(O)
  w <- colSums(O)
  E <- m %o% w
  Z <- O - E
  return(list(Z = Z, O = O, E = E, m = m, w = w))

}

#' @export
#'
#' @title Make data disjunctive
#'
#' @description Transforms a data matrix into a complete disjunctive matrix containing zeros and ones.
#'
#' @param DATA a matrix with the assumption that the columns contain categorical data.
#' @param impute_NA_to_mean a boolean flag on whether to impute missing values to their column means. Default: TRUE.
#'
#' @return
#' \item{DATAOUT}{Matrix: the expanded complete disjunctive data matrix that represents every level of every variable.}
#'
#'

make_data_disjunctive <- function (DATA, impute_NA_to_mean = TRUE)
{
  if (is.null(colnames(DATA))) {
    warning("'colnames(DATA)' were NULL. Setting to 1:ncol(DATA).")
    colnames(DATA) <- as.character(1:ncol(DATA))
  }
  data_dims <- dim(DATA)
  var_names <- colnames(DATA)
  ind_names <- rownames(DATA)
  new.col.count <- sum(apply(DATA, 2, function(x) {
    uniq <- unique(x)
    length(uniq) - sum(is.na(uniq))
  }))
  DATAOUT <- matrix(0, nrow(DATA), new.col.count)
  beginner <- 0
  variable.map <- new_colnames <- matrix(0, 1, 0)
  for (i in 1:data_dims[2]) {
    unique_elements <- unique(DATA[, i])
    unique_no_na <- unique_elements[!is.na(unique_elements)]
    mini.mat <- matrix(0, data_dims[1], length(unique_no_na))
    for (j in 1:ncol(mini.mat)) {
      mini.mat[which(DATA[, i] == unique_no_na[j]),
               j] <- 1
      new_colnames <- cbind(new_colnames, paste(var_names[i],
                                                ".", unique_no_na[j], sep = ""))
      variable.map <- cbind(variable.map, var_names[i])
    }
    these.missing <- which(rowSums(mini.mat) == 0)
    if (length(these.missing) > 0) {
      if (impute_NA_to_mean) {
        replacement <- colSums(mini.mat)/sum(colSums(mini.mat))
      }
      else {
        replacement <- rep(NA, ncol(mini.mat))
      }
      mini.mat[these.missing, ] <- matrix(replacement,
                                          length(these.missing), length(replacement),
                                          byrow = T)
    }
    ender <- beginner + length(unique_no_na)
    DATAOUT[, (beginner + 1):ender] <- mini.mat
    beginner <- ender
  }
  colnames(DATAOUT) <- c(new_colnames)
  rownames(DATAOUT) <- ind_names
  DATAOUT <- as.matrix(DATAOUT)
  attributes(DATAOUT)$variable.map <- c(variable.map)
  return(DATAOUT)
}

#' @export
#'
#' @title "Escofier's coding": a data doubling approach for continuous data
#'
#' @description A procedure to transform ("double") assumed continuous data into a format for use with Correspondence Analysis
#'
#' @param DATA A matrix of data with assumed continuous data for all columns
#' @param center A parameter to pass through to the \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{DATA}.
#' @param scale A parameter to pass through to the \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{DATA}.
#'
#' @return
#' \item{DATAOUT}{Matrix: of size \code{nrow(DATA)} rows and \code{ncol(DATA)*2} columns.}
#'

escofier_coding <- function (DATA, center = T, scale = T)
{
  if (is.null(colnames(DATA))) {
    warning("'colnames(DATA)' were NULL. Setting to 1:ncol(DATA).")
    colnames(DATA) <- as.character(1:ncol(DATA))
  }
  DATAOUT <- scale(DATA, center = center, scale = scale)
  dat.col.names <- c(paste0(colnames(DATAOUT), "-"), paste0(colnames(DATAOUT), "+"))
  DATA <- cbind((1 - DATAOUT)/2, (1 + DATAOUT)/2)
  colnames(DATAOUT) <- dat.col.names
  DATAOUT <- as.matrix(DATAOUT)
  attributes(DATAOUT)$variable.map <- gsub("\\-", "", gsub("\\+", "", dat.col.names))
  return(DATAOUT)
}


#' @export
#'
#' @title "Thermometer coding": a data doubling approach for ordinal data
#'
#' @description DESCRIBE HERE THAT FOR NOW IT ASSUMES THE ACTUAL MIN AND MAX ARE AVAILABLE IN THE COLUMNS.
#'
#' @param DATA A matrix of data with assumed ordinal data for all columns
#'
#' @return
#' \item{DATAOUT}{Matrix: of size \code{nrow(DATA)} rows and \code{ncol(DATA)*2} columns.}
#'

thermometer_coding <- function (DATA)
{
  if (is.null(colnames(DATA))) {
    warning("'colnames(DATA)' were NULL. Setting to 1:ncol(DATA).")
    colnames(DATA) <- as.character(1:ncol(DATA))
  }

  DATAOUT <- apply(DATA, 2, function(x) {
    max(x, na.rm = T) - x
  })

  mins <- apply(DATAOUT, 2, min, na.rm = T)
  maxs <- apply(DATAOUT, 2, max, na.rm = T)

  dat.col.names <- c(paste0(colnames(DATAOUT), "+"), paste0(colnames(DATAOUT), "-"))

  DATAOUT <- as.matrix(cbind(sweep(sweep(DATAOUT, 2, maxs, "-") * -1, 2, maxs, "/"), sweep(sweep(DATAOUT, 2, mins, "-"), 2, maxs, "/")))

  colnames(DATAOUT) <- dat.col.names
  attributes(DATAOUT)$variable.map <- gsub("\\-", "", gsub("\\+", "", dat.col.names))
  return(DATAOUT)
}


# this will take in the results of any of the methods and rebuild the X & Y (when possible) for a given rank
# rebuild_something <- function(){
#
# }
