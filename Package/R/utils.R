# a big set of utility functions


ca_preproc <- function(DATA){
  
  O <- DATA/sum(DATA)
  m <- rowSums(O)
  w <- colSums(O)
  E <- m %o% w
  Z <- O - E
  return(list(m = m, w = w, Z = Z, O = O, E = E))
  
}


make_data_disjunctive <- function (datain, impute.NA.to.mean = T) 
{
  if (is.null(colnames(datain))) {
    warning("'colnames(datain)' were NULL. Setting to 1:ncol(datain).")
    colnames(datain) <- as.character(1:ncol(datain))
  }
  data_dims <- dim(datain)
  var_names <- colnames(datain)
  ind_names <- rownames(datain)
  new.col.count <- sum(apply(datain, 2, function(x) {
    uniq <- unique(x)
    length(uniq) - sum(is.na(uniq))
  }))
  dataout <- matrix(0, nrow(datain), new.col.count)
  beginner <- 0
  variable.map <- new_colnames <- matrix(0, 1, 0)
  for (i in 1:data_dims[2]) {
    unique_elements <- unique(datain[, i])
    unique_no_na <- unique_elements[!is.na(unique_elements)]
    mini.mat <- matrix(0, data_dims[1], length(unique_no_na))
    for (j in 1:ncol(mini.mat)) {
      mini.mat[which(datain[, i] == unique_no_na[j]), 
               j] <- 1
      new_colnames <- cbind(new_colnames, paste(var_names[i], 
                                                ".", unique_no_na[j], sep = ""))
      variable.map <- cbind(variable.map, var_names[i])
    }
    these.missing <- which(rowSums(mini.mat) == 0)
    if (length(these.missing) > 0) {
      if (impute.NA.to.mean) {
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
    dataout[, (beginner + 1):ender] <- mini.mat
    beginner <- ender
  }
  colnames(dataout) <- c(new_colnames)
  rownames(dataout) <- ind_names
  dataout <- as.matrix(dataout)
  attributes(dataout)$variable.map <- c(variable.map)
  return(dataout)
}


escofier_coding <- function (DATA, center = T, scale = "SS1") 
{
  if (is.null(colnames(DATA))) {
    warning("'colnames(DATA)' were NULL. Setting to 1:ncol(DATA).")
    colnames(DATA) <- as.character(1:ncol(DATA))
  }
  DATA <- expo.scale(DATA, center = center, scale = scale)
  dat.col.names <- c(paste0(colnames(DATA), "-"), paste0(colnames(DATA), "+"))
  DATA <- cbind((1 - DATA)/2, (1 + DATA)/2)
  colnames(DATA) <- dat.col.names
  DATA <- as.matrix(DATA)
  attributes(DATA)$variable.map <- gsub("\\-", "", gsub("\\+", "", dat.col.names))
  return(DATA)
}


## this could incorporate alternate max's but now is not the time.
thermometer_coding <- function (DATA) 
{
  if (is.null(colnames(DATA))) {
    warning("'colnames(DATA)' were NULL. Setting to 1:ncol(DATA).")
    colnames(DATA) <- as.character(1:ncol(DATA))
  }
  
  DATA <- apply(DATA, 2, function(x) {
    max(x, na.rm = T) - x
  })
  
  mins <- apply(DATA, 2, min, na.rm = T)
  maxs <- apply(DATA, 2, max, na.rm = T)
  
  dat.col.names <- c(paste0(colnames(DATA), "+"), paste0(colnames(DATA), "-"))
  
  DATA <- as.matrix(cbind(sweep(sweep(DATA, 2, maxs, "-") * -1, 2, maxs, "/"), sweep(sweep(DATA, 2, mins, "-"), 2, maxs, "/")))
  
  colnames(DATA) <- dat.col.names
  attributes(DATA)$variable.map <- gsub("\\-", "", gsub("\\+", "", dat.col.names))
  return(DATA)
}


