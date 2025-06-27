#' Update Cholesky factors
update_L <- function(L, Lambda, Z, init=TRUE){
  if ( init ){
    LLt <- Matrix::t(Z %*% Lambda) %*% (Z %*% Lambda) + diag(dim(Z)[2])
    return( Cholesky(LLt, LDL=FALSE) )
  } else {
    return( Matrix::update(L, Matrix::t(Z%*%Lambda), mult=1) )
  }
}
