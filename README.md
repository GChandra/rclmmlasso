# rclmmlasso
Regression-Calibrated LMM Lasso for linear mixed-effect models with error-prone covariates.

## Installation

You can install the latest version of `rclmmlasso` from github with:

``` r
devtools::install_github("GChandra/rclmmlasso")
```

## Usage

The `rclmmlasso` function can be used to implement the Regression-Calibrated (RC) LMM Lasso method. The `ic.rclmmlasso` function implements the RC Lasso with AIC/BIC-based tuning parameter selection. An example for the implementation of `ic.rclmmlasso` is provided in its documentation.
