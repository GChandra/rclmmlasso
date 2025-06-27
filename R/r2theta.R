#' Residual function
r2theta <- function(y, X, Z, Lambda, u, beta){
  return (
    norm(y - X%*%beta - Z%*%Lambda%*%u, "2")^2 + norm(u, "2")^2
  )
}
