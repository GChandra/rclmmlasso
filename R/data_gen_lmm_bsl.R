#' Generate Data from Linear Mixed-Effects Model with Constant Covariates.
#' 
#' @description
#' Generates covariates, as well as error-prone versions, and produces a response
#' vector based on a specified Linear Mixed-Effects Model. The generated data
#' and metadata are returned in a list. Covariates are assumed to be constant.
#' @param n Number of clusters (e.g. individuals, clinics, etc.).
#' @param m Number of repeated measures per cluster (assumed equal across clusters).
#' @param p1 Number of error-prone covariates.
#' @param p2 Number of error-free covariates.
#' @param q Number of variance components.
#' @param k Scalar constant or n-dimensional vector of the number of replicates
#' for each cluster. If a scalar is provided, the number of replicates is assumed
#' to be constant across clusters.
#' @param mu Population mean for the LMM
#' @param beta_t Scalar coefficient for the time effect.
#' @param beta1 p1-dimensional coefficient vector for error-prone covariates.
#' @param beta2 p2-dimensional coefficient vector for error-free covariates.
#' @param mu_x (p1+p2)-dimensional vector for the mean of the covariate distribution.
#' @param Sigma_x (p1+p2)x(p1+p2)-dimensional covariance matrix for covariates.
#' @param Omega_0 q x q dimensional covariance matrix for random effects.
#' @param Sigma_d1 p1 x p1 dimensional covariance matrix for measurement errors.
#' @param z Index vector of covariates to include in random effects matrix. -1
#' indicates a random intercept, and 0 indicates a random slope associated with
#' time. All other positive integers correspond to columns of the covariate matrix
#' X.
#' @param Sigma_e m x m dimensional covariance matrix for residual errors.
#' @param theta d-dimensional vector of variance components for the random effects
#' covariance matrix.
#' @return a list with generated data and metadata. Importantly, Xs is the generated
#' error-prone data of dimension sum(k) x (p1+p2). The first k_1 rows correspond to
#' the k_1 replicates for cluster 1, the next k_2 to the k_2 replicates for cluster
#' 2, and so on.
#' 
#' @export
#' 
data_gen_lmm_bsl <- function(n, m, p1, p2, q, k, mu, beta_t, beta1, beta2, mu_x,
                             Sigma_x, Omega_0, Sigma_d1, z, Sigma_e, theta){
  if (length(k)==1){
    k <- rep(k, n)
  }
  
  t <- rep(0:(m-1), n)
  
  X <- MASS::mvrnorm(n, mu_x, Sigma_x)
  
  Xs <- X[rep(1:nrow(X), times=k),]
  
  Delta <- cbind(
    MASS::mvrnorm(sum(k), rep(0, p1), Sigma_d1),
    matrix(0, nrow=sum(k), ncol=p2)
  )
  
  Xs <- Xs + Delta
  
  X <- X[rep(1:nrow(X), each=m),]
  
  Sigma_d <- rbind( cbind(Sigma_d1, matrix(0, p1, p2)), matrix(0, p2, p1+p2) )
  
  b <- MASS::mvrnorm(n, rep(0, q), Omega_0)
  b <- c(t(b))
  
  z_ord <- order(z)
  z <- sort(z)
  Z <- list()
  for (i in 1:n){
    inds <- ((i-1)*m + 1) : (i*m)
    
    # intercept
    if (z[1] == -1) {
      # intercept + time
      if (length(z)>1 && z[2] == 0) {
        Z[[i]] <- cbind(
          rep(1,m),
          t[inds],
          X[inds, z[-(1:2)]]
        )
        # intercept, no time
      } else { 
        Z[[i]] <- cbind(
          rep(1,m),
          X[inds, z[-1]]
        )
      }
      # time, no intercept
    } else if (z[1] == 0){
      Z[[i]] <- cbind(
        t[inds],
        X[inds, z[-1]]
      )
      # no intercept, no time
    } else {
      Z[[i]] <- X[inds, z]
    }
    
    Z[[i]] <- Z[[i]][, z_ord]
  }
  Z <- bdiag(Z)
  
  e <- MASS::mvrnorm(n, rep(0, m), Sigma_e)
  e <- c(t(e))
  
  beta_x <- c(beta1, beta2)
  y <- c(as.matrix(
    mu + t*beta_t + X %*% beta_x + Z %*% b + e
  ))
  
  beta <- c(beta_t, beta_x)
  S <- which(beta != 0)
  
  return (
    list(t=t, X=X, Xs=Xs, Delta=Delta, Z=Z, b=b, e=e, y=y,
         cluster=as.factor(rep(1:n, each=m)), beta_t=beta_t, beta_x=beta_x,
         beta1=beta1, beta2=beta2, beta=beta, S=S, mu=mu,
         mu_x=mu_x, Sigma_x=Sigma_x,
         Omega_0=Omega_0, Sigma_d1=Sigma_d1, Sigma_d=Sigma_d,
         Sigma_e=Sigma_e, theta=theta, n=n, m=m, N=n*m, p1=p1, p2=p2, p_x=p1+p2,
         p=1+p1+p2, q=q, k=k, d=length(theta), n_obs=n*m,
         n_params=2+p1+p2+1+length(theta))
  )
}
