<!-- badges: start -->
  [![Travis build status](https://travis-ci.com/lasy/HiddenSemiMarkov.svg?branch=master)](https://travis-ci.com/lasy/HiddenSemiMarkov)
  <!-- badges: end -->
  
  :warning: The package is not operational yet. It will be up and running shortly :warning:

# HiddenSemiMarkov

`HiddenSemiMarkov` is an R package for specifying hidden semi-Markov models, simulating sequences, fitting models to specific observation sequences and predicting the sequence of hidden states from observation sequences.

It is especially designed to model multivariate sequences with frequent missing data. The probabilities of missing data may be state- or variable-dependent and specified manually or learned from observation sequences.

## Installation

To install the package from github:
* Open R or Rstudio
* Make sure that the `devtools` library is installed or install it with:

`install.packages("devtools")`

* Load the "devtools" library with: 

`library(devtools)`

* Install the `HiddenSemiMarkov` package with the command:

`devtools::install_github("lasy/HiddenSemiMarkov")` 

or, if you wish to build the vignette for this package:

`devtools::install_github("lasy/HiddenSemiMarkov", build_vignettes = TRUE)` 


## To get started

Find examples in the vignette by typing `browseVignettes(package = "HiddenSemiMarkov")` 
or start specifying a model with the `specify_hsmm` function (type `?specify_hsmm` for documentation).
