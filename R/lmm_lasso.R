#' Coordinate-descent LMM Lasso Algorithm
#' 
#' @description Coordinate-descent based algorithm to fit the LMM Lasso.
#' @param y response vector
#' @param X covariate matrix
#' @param Z random effects matrix
#' @param theta_0 initial estimate for variance components vector (optional)
#' @param beta_0 initial estimate for fixed-effects coefficient vector (optional)
#' @param sig2e_0 initial estimate for residual error variance (optional)
#' @param p dimension of covariate matrix X (also, number of fixed effects)
#' @param q number of random effects
#' @param d variance components vector dimension
#' @param n number of clusters
#' @param N total number of response observations, including repeated measures
#' @param REML specify whether the ML or REML estimator will be used (default TRUE)
#' @param BIC specify whether BIC or AIC will be returned for the fitted model.
#' @param lambda scalar penalty parameter
#' @param penalty.factor Optional p-dimensional vector of penalty weights, so that
#' different predictors can have different penalties. Default is equal weighting, and
#' weights are normalized to sum to 1.
#' @param eta Threshold for convergence of outer loop (which solves the mixed model equations)
#' @param eta_gs Threshold for convergence of inner loop (which performs coordinate descent)
#' @param maxits Maximum iterations before exit
#' @param standardize specify whether to standardize predictors (i.e. mean-center and
#' scale to have unit variance) for the fitting procedure. Final coefficient estimates
#' will be rescaled to the original predictor scales.
#' @param intercept specify whether the model should include an intercept.
#' @return list with eta (convergence value), AIC/BIC value, variance component
#' estimate, random effects covariance matrix estimate, intercept, fixed-effects,
#' and random-effects estimates, residual error variance estimate, and total
#' number of iterations for the optimal model.
#' @export
lmm_lasso <- function(y, X, Z,
                      theta_0=NULL, beta_0=NULL, sig2e_0=NULL,
                      p, q, d, n, N,
                      REML=TRUE, BIC=TRUE, lambda, penalty.factor=NULL,
                      eta=1e-4, eta_gs=1e-3, maxits=1000,
                      standardize=TRUE, intercept=TRUE){
  # center data
  y0 <- as.matrix(y) - mean(y)*intercept
  X0 <- scale(X, center=intercept, scale=standardize)
  
  if (is.null(penalty.factor)){
    penalty.factor <- rep(1, p)
  } else {
    # rescale factors to add to p, as in glmnet
    penalty.factor <- penalty.factor * p / sum(penalty.factor)
  }
  
  if (is.null(theta_0) | is.null(beta_0) | is.null(sig2e_0)){
    beta_k <- rep(0, p)
    theta_k <- rep(1, d)
    sig2e_k <- 1
  } else {
    beta_k <- beta_0
    theta_k <- theta_0
    sig2e_k <- sig2e_0
  }
  
  # create lower bound vector for nloptr (lower-triangular reference is based on
  # the choice in the mapping() function).
  lb <- matrix(-Inf, nrow=q, ncol=q)
  diag(lb) <- 0
  lb <- lb[lower.tri(lb, diag=TRUE)]
  
  # create covariance Cholesky factor
  Lambda_k <- mapping(theta_k, n, q, init=TRUE)
  L <- update_L(L=NULL, Lambda=Lambda_k, Z=Z, init=TRUE)
  
  u_k <- rep(0, n*q)
  d_k <- c(u_k, beta_k)
  
  d_old <- d_k + 2*eta
  
  k <- 0
  while ( (max(abs(d_k-d_old)) > eta) & (k < maxits) ) {
    d_old <- d_k
    
    # update covariance Cholesky factor
    Lambda_k@x[] <- mapping(theta_k, n, q, init=FALSE)
    L <- update_L(L=L, Lambda=Lambda_k, Z=Z, init=FALSE)
    
    # Solve normal equations
    R_ZX <- Matrix::solve( L, Matrix::crossprod(Z%*%Lambda_k, X0), system="L" )
    R_XtR_X <- as( Matrix::crossprod(X0) - Matrix::crossprod(R_ZX), "dpoMatrix" )
    
    # outer loop for convergence
    # set active-set
    J <- which(beta_k != 0)
    if ( (k < 20) || (k %% 200 == 0) ){
      J <- 1:p
    }
    
    d_k_i <- c(u_k, beta_k)
    d_old_i <- d_k_i + 2*eta_gs
    while(max(abs(d_k_i-d_old_i)) > eta_gs){
      d_old_i <- d_k_i
      
      u_k <- c(as.matrix(
        Matrix::solve(L, Matrix::crossprod(Z%*%Lambda_k, y0) -
                        Matrix::crossprod(Z%*%Lambda_k, X0%*%beta_k),
              system="LDLt")
      ))
      
      max_diff <- 2*eta_gs
      while (max_diff > eta_gs){
        max_diff <- -Inf
        
        for (j in J){
          beta_j_new <- st( t(X0[,j]) %*%
                              (y0 - X0[,-j]%*%beta_k[-j] - Z%*%Lambda_k%*%u_k) /
                              norm(X0[,j], "2")^2,
                            lambda * penalty.factor[j] / norm(X0[,j], "2")^2)
          
          if (abs(beta_j_new - beta_k[j]) > max_diff){
            max_diff <- abs(beta_j_new - beta_k[j])
          }
          
          beta_k[j] <- beta_j_new
        }
      }
      
      d_k_i <- c(u_k, beta_k)
    }
    
    d_k <- d_k_i
    
    # Update residual
    if (REML){
      sig2e_k <- r2theta(y0, X0, Z, Lambda_k, u_k, beta_k) / (N-p)
    } else {
      sig2e_k <- r2theta(y0, X0, Z, Lambda_k, u_k, beta_k) / N
    }
    
    # optimization step
    opt_dev <- nloptr(x0=theta_k, eval_f=devtheta, eval_grad_f=NULL,
                      lb=lb,
                      opts=list("algorithm"="NLOPT_LN_BOBYQA",
                                "xtol_rel"=1.0e-4),
                      y=y0, X=X0, Z=Z, L=L, Lambda=Lambda_k,
                      sig2e=sig2e_k, uk=u_k,
                      beta=beta_k,
                      n=n, N=N, p=p, q=q,
                      R_XtR_X=R_XtR_X, REML=REML)
    # adding 1e-10 to avoid issues with zeroing out elements I need to access.
    theta_k <- opt_dev$solution
    theta_k <- theta_k + 1e-10*sign(theta_k)
    
    k <- k + 1
  }
  
  if (k == maxits){ warning("algorithm failed to converge, reached max iterations.") }
  
  # update covariance Cholesky factor and u_k once more
  Lambda_k@x[] <- mapping(theta_k, n, q, init=FALSE)
  L <- update_L(L=L, Lambda=Lambda_k, Z=Z, init=FALSE)
  # Solve normal equations
  R_ZX <- Matrix::solve( L, Matrix::crossprod(Z%*%Lambda_k, X0), system="L" )
  R_XtR_X <- as( Matrix::crossprod(X0) - Matrix::crossprod(R_ZX), "dpoMatrix" )
  u_k <- c(as.matrix(
    Matrix::solve(L, Matrix::crossprod(Z%*%Lambda_k, y0) -
                    Matrix::crossprod(Z%*%Lambda_k, X0%*%beta_k),
          system="LDLt")
  ))
  
  if (standardize){
    beta_k <- beta_k / attr(X0, "scaled:scale")
  }
  
  if (BIC) {
    ic <- opt_dev$objective + log(N)*(sum(beta_k != 0) + d)
  } else {
    ic <- opt_dev$objective + 2*(sum(beta_k != 0) + d)
  }
  
  mu_k <- 0
  if (intercept){
    mu_k <- mean(y - X %*% beta_k - Z %*% Lambda_k %*% u_k)
  }
  Omega_0_opt <- as.matrix(Lambda_k[1:q,1:q]) %*% t(as.matrix(Lambda_k[1:q,1:q])) * sig2e_k
  
  return (list(
    eta_opt=max(abs(d_k-d_old)),
    ic_opt=ic,
    theta_opt=theta_k,
    Omega_0_opt=Omega_0_opt,
    mu_opt=mu_k,
    beta_opt=beta_k,
    u_opt=u_k,
    sig2e_opt=sig2e_k,
    iter=k
  ))
}
