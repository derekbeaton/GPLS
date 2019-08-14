# a big set of utility functions


ca_preproc <- function(DATA){
  
  O <- DATA/sum(DATA)
  m <- rowSums(O)
  w <- colSums(O)
  E <- m %o% w
  Z <- O - E
  return(list(m = m, w = w, Z = Z, O = O, E = E))
}

