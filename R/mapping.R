#' Mapping of variance components to random effects covariance matrix
mapping <- function(theta, n, q, init=TRUE){
  if ( init ){
    Lambda_0 <- matrix(0, q, q)
    Lambda_0[lower.tri(Lambda_0, diag=TRUE)] <- theta
    Lambda <- Matrix::bdiag( rep(list(Lambda_0), n) )
    
    return( Lambda )
  } else {
    return ( rep(theta, n) )
  }
}
