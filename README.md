---
editor_options: 
  markdown: 
    wrap: 72
---

# Dirichlet Process Mixture Models in R

This repository demonstrates the application of Bayesian nonparametric
methods using the `dirichletprocess` package in R. The Dirichlet process
mixture model (DPMM) is a framework for density estimation and
clustering that automatically determines the appropriate number of
mixture components from the data.

## Overview

The `dirichletprocess` package provides a flexible framework for fitting
complex Bayesian nonparametric models in R. Unlike parametric models
that require a priori specification of the number of components, DPMMs
determine this automatically, making them ideal for exploratory data
analysis and when the true number of clusters is unknown.

## Tasks

This repository includes three tasks.

### Easy: Standard Mixture Models

Application of built-in mixture models to standard datasets: - Fitting a
Gaussian mixture to the `faithful` dataset (univariate) - Fitting a
multivariate normal mixture to the `iris` dataset - Visualization and
interpretation of the results

### Medium: Lognormal Mixture Analysis

Exploration of how the alpha parameter affects clustering behavior: -
Generation of synthetic data from a lognormal mixture - Fitting models
with different alpha priors - Analysis of posterior samples and
convergence diagnostics - Visualization of cluster numbers and their
relationship to prior choice

### Hard: Custom Gamma Mixture Model

Implementation of a custom mixture model not included in the package: -
Development of a Gamma mixture model with conjugate sampling -
Mathematical framework and implementation details - Testing on simulated
data - Comparison with true generating distribution

## Usage

The main analysis is contained in the R Markdown document
`dirichletprocess_analysis.Rmd`.

``` r
rmarkdown::render("dirichletprocess_analysis.Rmd")
```

Individual R scripts for each task are: - `easy.R` - `medium.R` -
`hard.R`

## Mathematical Framework

The Dirichlet process mixture model can be represented as:

$$p(x) = \sum_{j=1}^{\infty} w_j \cdot k(x | \theta_j)$$

where: - $w_j$ are the mixture weights - $k(x | \theta_j)$ is the kernel
density with parameters $\theta_j$ - $\theta_j \sim G$ where
$G \sim DP(\alpha, G_0)$

The concentration parameter $\alpha$ controls the expected number of
clusters, while the base measure $G_0$ represents our prior belief about
the parameters.

## References

1.  Ross, G. J., & Markwick, D. (2018). dirichletprocess: An R Package
    for Fitting Complex Bayesian Nonparametric Models.
    [arXiv:1809.07734](https://arxiv.org/abs/1809.07734)

2.  Neal, R. M. (2000). Markov Chain Sampling Methods for Dirichlet
    Process Mixture Models. Journal of Computational and Graphical
    Statistics, 9(2), 249-265.

3.  Escobar, M. D., & West, M. (1995). Bayesian Density Estimation and
    Inference Using Mixtures. Journal of the American Statistical
    Association, 90(430), 577-588.
