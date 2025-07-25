% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lmm_lasso.R
\name{lmm_lasso}
\alias{lmm_lasso}
\title{Coordinate-descent LMM Lasso Algorithm}
\usage{
lmm_lasso(
  y,
  X,
  Z,
  theta_0 = NULL,
  beta_0 = NULL,
  sig2e_0 = NULL,
  p,
  q,
  d,
  n,
  N,
  REML = TRUE,
  BIC = TRUE,
  lambda,
  penalty.factor = NULL,
  eta = 1e-04,
  eta_gs = 0.001,
  maxits = 1000,
  standardize = TRUE,
  intercept = TRUE
)
}
\arguments{
\item{y}{response vector}

\item{X}{covariate matrix}

\item{Z}{random effects matrix}

\item{theta_0}{initial estimate for variance components vector (optional)}

\item{beta_0}{initial estimate for fixed-effects coefficient vector (optional)}

\item{sig2e_0}{initial estimate for residual error variance (optional)}

\item{p}{dimension of covariate matrix X (also, number of fixed effects)}

\item{q}{number of random effects}

\item{d}{variance components vector dimension}

\item{n}{number of clusters}

\item{N}{total number of response observations, including repeated measures}

\item{REML}{specify whether the ML or REML estimator will be used (default TRUE)}

\item{BIC}{specify whether BIC or AIC will be returned for the fitted model.}

\item{lambda}{scalar penalty parameter}

\item{penalty.factor}{Optional p-dimensional vector of penalty weights, so that
different predictors can have different penalties. Default is equal weighting, and
weights are normalized to sum to 1.}

\item{eta}{Threshold for convergence of outer loop (which solves the mixed model equations)}

\item{eta_gs}{Threshold for convergence of inner loop (which performs coordinate descent)}

\item{maxits}{Maximum iterations before exit}

\item{standardize}{specify whether to standardize predictors (i.e. mean-center and
scale to have unit variance) for the fitting procedure. Final coefficient estimates
will be rescaled to the original predictor scales.}

\item{intercept}{specify whether the model should include an intercept.}
}
\value{
list with eta (convergence value), AIC/BIC value, variance component
estimate, random effects covariance matrix estimate, intercept, fixed-effects,
and random-effects estimates, residual error variance estimate, and total
number of iterations for the optimal model.
}
\description{
Coordinate-descent based algorithm to fit the LMM Lasso.
}
