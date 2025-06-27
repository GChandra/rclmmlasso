#' Choose Covariance Matrix Structure
#'
#' @description Returns a covariance matrix of the specified dimension and structure.
#' @param sig2 Scalar variance of each covariate (assumed constant across covariates).
#' @param rho Correlation for AR-1 structure (and block-diagonal structure with AR-1
#' blocks).
#' @param p Covariance matrix dimension.
#' @param blk_sz Size of blocks for block-diagonal structure.
#' @param structure Covariance matrix structure (diag, block, or AR1).
#' @return Covariance matrix
#' 
#' @export
#' 
choose_matrix <- function(sig2, rho=0.7, p, blk_sz = 20, structure){
  if(structure == "diag") {
    mat <- matrix(0, p, p)
    diag(mat) <- sig2
  } else if(structure == "block") {
    blk <- matrix(NA, blk_sz, blk_sz)
    for (i in 1:blk_sz){
      for (j in 1:blk_sz){
        blk[i,j] <- sig2 * rho^(abs(i-j))
      }
    }
    mat <- as.matrix( Matrix::bdiag(replicate(p/blk_sz, blk, simplify=FALSE)) )
  } else if(structure=="AR1"){
    mat <- matrix(0, p, p)
    for (i in 1:p){
      for (j in 1:p){
        mat[i,j] <- sig2 * rho^abs(i-j)
      }
    }
  } else {
    mat <- matrix(NA, p, p)
  }
  
  return(mat)
}
