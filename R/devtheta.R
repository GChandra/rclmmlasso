#' Deviance function (twice negative log-likelihood) parametrized by variance
#' components.
devtheta <- function(theta, y, X, Z, L, Lambda, sig2e, uk, beta, n, N, p, q,
                     R_XtR_X=NULL, REML=TRUE){
  u <- uk
  if (REML & is.null(R_XtR_X)){
    REML <- FALSE
    message("R_X matrix not provided; cannot do REML optimization.
            Switching to normal optimization.")
  }
  
  Lambda@x[] <- mapping(theta, n, q, init=FALSE)
  L <- update_L(L, Lambda, Z, init=FALSE)
  
  r2theta <- r2theta(y, X, Z, Lambda, u, beta)
  
  if (REML){
    return(
      2*Matrix::determinant(L, sqrt=TRUE, logarithm=TRUE)$modulus +
        Matrix::determinant(R_XtR_X, logarithm=TRUE)$modulus +
        (N-p)*(1 + log(2*pi*r2theta/(N-p)))
    )
  } else {
    return(
      2*Matrix::determinant(L, sqrt=TRUE, logarithm=TRUE)$modulus + N*(1 + log(2*pi*r2theta/N))
    )
  }
}
