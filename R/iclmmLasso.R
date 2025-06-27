#' Coordinate-descent LMM Lasso Algorithm with AIC/BIC-based Tuning Parameter Selection
#' 
#' @description Coordinate-descent based algorithm to fit the LMM Lasso with AIC/BIC-based tuning parameter selection.
#' @param formula LMM formula to fit
#' @param data data frame which includes the variables from the formula
#' @param p dimension of covariate matrix X (also, number of fixed effects)
#' @param q number of random effects
#' @param d variance components vector dimension
#' @param n number of clusters
#' @param N total number of response observations, including repeated measures
#' @param theta_0 initial estimate for variance components vector (optional)
#' @param beta_0 initial estimate for fixed-effects coefficient vector (optional)
#' @param sig2e_0 initial estimate for residual error variance (optional)
#' @param lambda.max maximum value for penalty parameter over which to optimize.
#' If not provided, one will be computed using a data-driven approach.
#' @param nlambda Number of penalty parameters to search through.
#' @param lambda.min.ratio Ratio which will determine minimum lambda (lambda min =
#' lambda.max*lambda.min.ratio). Default is 0.001.
#' @param penalty.factor Optional p-dimensional vector of penalty weights, so that
#' different predictors can have different penalties. Default is equal weighting, and
#' weights are normalized to sum to 1.
#' @param REML specify whether the ML or REML estimator will be used (default TRUE)
#' @param BIC specify whether BIC or AIC will be returned for the fitted model.
#' @param eta Threshold for convergence of outer loop (which solves the mixed model equations)
#' @param eta_gs Threshold for convergence of inner loop (which performs coordinate descent)
#' @param maxits Maximum iterations before exit
#' @param standardize specify whether to standardize predictors (i.e. mean-center and
#' scale to have unit variance) for the fitting procedure. Final coefficient estimates
#' will be rescaled to the original predictor scales.
#' @param progress 0, 1, or 2. Larger integers will result in printing more information
#' on the progress of the algorithm.
#' @param getInit specify whether you would like just the data driven choice of lambda
#' and the computed initial values for theta, beta, and sig2e to be returned. If so,
#' the algorithm will not be run to completion; only the initial values will be
#' computed and returned. Default is FALSE.
#' @param ... extra parameters for lmer function
#' @return list with AIC/BIC value, variance component estimate, random effects
#' covariance matrix estimate, intercept, fixed-effects, and random-effects
#' estimates, residual error variance estimate, optimal penalty parameter, and
#' matrix of indicators S_hat. S_hat is an nlambda x p matrix with 1's where the
#' model produced a nonzero coefficient estimate, and 0's elsewhere.
#' @export
ic.lmm_lasso <- function(formula, data,
                          p, q, d, n, N,
                          theta_0=NULL, beta_0=NULL, sig2e_0=NULL,
                          lambda.max=NULL, nlambda=20, lambda.min.ratio=0.001,
                          penalty.factor=NULL, REML=TRUE, BIC=TRUE,
                          eta=1e-4, eta_gs=5e-3, maxits=5e4,
                          standardize=TRUE, progress=FALSE,
                          getInit=FALSE, ...){
  if (is.null(penalty.factor)){
    penalty.factor <- rep(1, p)
  } else {
    # rescale factors to sum to p, as in glmnet
    penalty.factor <- penalty.factor * p / sum(penalty.factor)
  }
  
  formula <- formula(formula)
  
  frm <- lFormula(eval(formula), data)
  Z <- Matrix::t(frm$reTrms$Zt)
  
  intercept <- all(frm$X[,1]==1)
  y <- frm$fr[,1]
  if (intercept){
    X <- frm$X[,-1]
  } else {
    X <- frm$X
  }
  
  y0 <- y - mean(y)*intercept
  X0 <- scale(X, center=intercept, scale=standardize)
  
  if (is.null(lambda.max)){
    if(any(penalty.factor==0)){
      fit <- lmer(
        paste0(frm$formula[[2]], " ~ ", as.numeric(intercept), " + ",
               paste(colnames(frm$X)[-1][penalty.factor==0], collapse=" + "), " + ",
               paste("(", findbars(formula), ")", collapse=" + ")
        ),
        data,
        ...
      )
    } else {
      fit <- lmer(
        paste0(frm$formula[[2]], " ~ ", as.numeric(intercept), " + ",
               paste("(", findbars(formula), ")", collapse=" + ")
        ),
        data,
        ...
      )
    }
    
    beta_0 <- rep(0, p)
    beta_0[penalty.factor==0] <- fixef(fit)[-1]
    sig2e_0 <- as.data.frame(VarCorr(fit))$vcov[d+1]
    theta_0 <- c(vech(t(chol(as.matrix(Matrix::bdiag(VarCorr(fit))))))) / sqrt(sig2e_0)
    
    lambda.max.ub <- max(abs( t(X0) %*% y0 )) * sum(penalty.factor) / p
    fit <- lmm_lasso(y, X, Z,
                     theta_0, beta_0, sig2e_0,
                     p, q, d, n, N,
                     REML=REML, BIC=BIC,
                     lambda=lambda.max.ub, penalty.factor=penalty.factor,
                     eta=eta, eta_gs=eta_gs, maxits=maxits,
                     standardize=standardize, intercept=intercept)
    theta_0 <- fit$theta_opt
    beta_0 <- fit$beta_opt
    sig2e_0 <- fit$sig2e_opt
    
    # create covariance Cholesky factor
    Lambda_0 <- mapping(theta_0, n, q)
    
    lambda.max <- max(abs( t(X0) %*% (y0 -
                                        X0%*%beta_0 - Z%*%Lambda_0%*%fit$u_opt) )) *
      sum(penalty.factor) / p
    
    if (getInit){
      return (list(
        lambda.max=lambda.max,
        theta_0=theta_0,
        beta_0=beta_0,
        sig2e_0=sig2e_0
      ))
    }
  }
  
  lambda.seq <- exp(
    seq( log(lambda.max), log(lambda.min.ratio*lambda.max), length.out=nlambda )
  )
  
  if (progress){
    start <- proc.time()
  }
  S_hat <- matrix(NA, nrow=nlambda, ncol=p)
  ic_min <- Inf
  for (lambda in lambda.seq){
    if (progress){
      print(lambda)
    }
    
    fit <- lmm_lasso(y, X, Z,
                     theta_0, beta_0, sig2e_0,
                     p, q, d, n, N,
                     REML=REML, BIC=BIC,
                     lambda=lambda, penalty.factor=penalty.factor,
                     eta=eta, eta_gs=eta_gs, maxits=maxits,
                     standardize=standardize, intercept=intercept)
    theta_0 <- fit$theta_opt
    beta_0 <- fit$beta_opt
    
    if (progress==2){
      beta_p <- beta_0
      names(beta_p) <- NULL
      print(beta_p)
      print(theta_0)
      print(fit$ic_opt)
    }
    if (standardize){
      beta_0 <- beta_0 * attr(X0, "scaled:scale")
    }
    
    sig2e_0 <- fit$sig2e_opt
    
    S_hat[ which(lambda.seq==lambda) , ] <- as.numeric(beta_0 != 0)
    
    if (fit$ic_opt < ic_min){
      ic_min <- fit$ic_opt
      
      theta_opt <- fit$theta_opt
      mu_opt <- fit$mu_opt
      beta_opt <- fit$beta_opt
      u_opt <- fit$u_opt
      sig2e_opt <- fit$sig2e_opt
      Omega_0_opt <- fit$Omega_0_opt
      lambda_opt <- lambda
    }
  }
  if (progress){
    end <- proc.time()
    print(paste( round((end-start)[3], 2), "seconds" ))
  }
  
  return (list(
    ic_min=ic_min,
    mu_opt=mu_opt,
    beta_opt=beta_opt,
    u_opt=u_opt,
    S_hat=S_hat,
    sig2e_opt=sig2e_opt,
    theta_opt=theta_opt,
    Omega_0_opt=Omega_0_opt,
    lambda_opt=lambda_opt
  ))
}
