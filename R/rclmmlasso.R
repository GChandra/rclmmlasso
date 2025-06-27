#' Regression-Calibrated LMM Lasso with AIC/BIC Tuning Parameter Selection
#' 
#' @description RC LMM Lasso performs Lasso for LMMs with RC correction for
#' error-prone data.
#' @param formula LMM formula with variables from the data provided.
#' @param data Data (list) required for the LMM. At minimum, this list should
#' include the following elements: Xs, the sum(k) x (p1+p2) dimensional covariate
#' matrix with error-prone data; n, the number of clusters; k, a vector or scalar
#' with the number of replicates per cluster; p1, the number of error-free covariates;
#' p2, the number of error-prone covariates; m, the number of repeated measures
#' (assumed equal across clusters); t, the time vector; y, the response vector;
#' cluster, an m*n-dimensional vector of cluster indicators. The names of the columns
#' of Xs should include the variables in the formula.
#' @param ... Remaining arguments for the CD LMM Lasso (ic.lmm_lasso) function.
#' @return list results from running the CD LMM Lasso method.
#' @examples
#' library(matrixcalc)
#' library(nloptr)
#' library(lme4)
#' library(tidyr)
#' library(dplyr)
#' 
#' set.seed(281)
#' 
#' n <- 1000
#' m <- 3
#' k <- sample(2:3, size=n, replace=TRUE)
#' p1 <- 30
#' p2 <- 2
#' q <- 1
#' s1 <- 5
#'
#' beta_t <- 0
#'
#' S1 <- sample(1:p1, size=s1)
#' beta1 <- rep(0, p1)
#' beta1[S1] <- rnorm(s1, mean=0.8, sd=1)
#' beta2 <- rnorm(2, mean=0.8, sd=1)
#' 
#' sig2e <- 0.6
#' Sigma_e <- sig2e * diag(m)
#' 
#' sig2b <- 15
#' Omega_0 <- matrix(c(sig2b),
#'                   nrow=q, ncol=q)
#' Omega_0_c <- t(chol(Omega_0)) / sqrt(sig2e)
#' theta <- Omega_0_c[lower.tri(Omega_0_c, diag=TRUE)]
#' 
#' z <- c(-1)
#' 
#' mu <- 20
#' mu_x <- c(rep(4, p1+p2))
#' mu_x[(p1+1):(p1+p2)] <- c(63, 0.5)
#' 
#' sig2x1 <- 1
#' Sigma_x1x1 <- choose_matrix(sig2=sig2x1, rho=0.3, p=p1, structure="AR1")
#' Sigma_x1x2 <- matrix(0, nrow=p1, ncol=p2)
#' Sigma_x1x2[1,2] <- -0.2*sqrt(0.25)
#' Sigma_x2x2 <- matrix(c(35, 0, 0, 0.25), nrow=p2, ncol=p2)
#' 
#' Sigma_x <- rbind(
#'   cbind(Sigma_x1x1, Sigma_x1x2),
#'   cbind(t(Sigma_x1x2), Sigma_x2x2)
#' )
#' 
#' sig2d1 <- sig2x1
#' Sigma_d1 <- sig2d1 * diag(p1)
#' 
#' data <- data_gen_lmm_bsl(n=n, m=m, p1=p1, p2=p2, q=q, k=k, mu=mu,
#'                          beta_t=beta_t, beta1=beta1, beta2=beta2,
#'                          mu_x=mu_x, Sigma_x=Sigma_x, Omega_0=Omega_0,
#'                          Sigma_d1=Sigma_d1, z=z, Sigma_e=Sigma_e,
#'                          theta=theta)
#' colnames(data$Xs) <- paste0("X", 1:(data$p_x))
#' colnames(data$X) <- paste0("X", 1:(data$p_x))
#' 
#' formula <- paste("y ~ 1 + t + ", paste( colnames(data$Xs),
#'                                         collapse=" + "),
#'                  " + (1 | cluster)",
#'                  sep="" )
#' 
#' 
#' fit.ic.rclmmlasso <- ic.rclmmlasso(formula, data=data,
#'                                    p=data$p, q=data$q, d=data$d,
#'                                    n=length(unique(data$cluster)), N=n*m,
#'                                    lambda.min.ratio=0.01,
#'                                    REML=TRUE, BIC=TRUE,
#'                                    standardize=TRUE, progress=2, getInit=FALSE,
#'                                    control=lmerControl(optCtrl=list(xtol_rel=1e-5)))
#' print("done")
#'                                    
#' @export
#' 
ic.rclmmlasso <- function(formula, data, ...){
  Sigma_u_hat <- data.frame(id=rep(1:data$n, times=data$k), data$Xs)
  Sigma_u_hat <- split( Sigma_u_hat , f = Sigma_u_hat$id )
  Sigma_u_hat <- lapply(
    Sigma_u_hat,
    FUN = function(X) t(scale(X[,-1], scale=FALSE)) %*% scale(X[,-1], scale=FALSE)
  )
  Sigma_u_hat <- Reduce("+", Sigma_u_hat) / sum(data$k-1)
  
  mu_x_hat <- mu_w_hat <- colMeans( data$Xs )
  
  nu <- sum(data$k) - (sum(data$k^2) / sum(data$k))
  
  Xs_bar <- data.frame(id=rep(1:data$n, times=data$k), data$Xs)
  Xs_bar <- Xs_bar %>%
    group_by(id) %>%
    summarise_all(mean)
  Xs_bar <- as.matrix(Xs_bar[,-1])
  
  Sigma_x_hat <- t(Xs_bar - rep(1, data$n) %*% t(mu_x_hat)) %*%
    diag(data$k) %*%
    (Xs_bar - rep(1, data$n) %*% t(mu_x_hat))
  Sigma_x_hat <- (Sigma_x_hat - (data$n-1)*Sigma_u_hat) / nu
  
  Xs_bar_1 <- data$Xs[ !duplicated(rep(1:data$n, times=data$k)), ]
  
  Lambda_hat <- solve(Sigma_x_hat+Sigma_u_hat, Sigma_x_hat)[,1:data$p1]
  X_hat <- rep(1,nrow(Xs_bar_1)) %*% (t(mu_w_hat[1:data$p1]) - ( mu_w_hat %*% Lambda_hat )) +
    Xs_bar_1 %*% Lambda_hat
  X_hat <- cbind(X_hat, Xs_bar_1[,(data$p1+1):(data$p1+data$p2)])
  
  dat <- data.frame(y=data$y,
                    t=data$t,
                    X_hat[rep(1:nrow(X_hat), each=data$m),],
                    cluster=data$cluster)
  
  fit.ic.rcu <- ic.lmm_lasso(formula, dat, ...)
  
  # de-bias variance components
  Omega_tau_0 <- matrix(c(
    t(fit.ic.rcu$beta_opt[-1]) %*% Sigma_x_hat %*%
      (diag(data$p_x) - solve(Sigma_x_hat+Sigma_u_hat, Sigma_x_hat)) %*%
      fit.ic.rcu$beta_opt[-1]
  ),
  nrow=data$q, ncol=data$q)
  
  Omega_0_RC_hat <- fit.ic.rcu$Omega_0_opt - Omega_tau_0
  
  fit.ic.rcu$Omega_0_RC_hat <- Omega_0_RC_hat
  
  return ( fit.ic.rcu )
}
