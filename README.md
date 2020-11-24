<!-- badges: start -->
  [![Travis build status](https://travis-ci.com/lasy/HiddenSemiMarkov.svg?branch=master)](https://travis-ci.com/lasy/HiddenSemiMarkov)
  <!-- badges: end -->

# HiddenSemiMarkov
HiddenSemiMarkov R package

To install the package from github:
* Open R or Rstudio
* Make sure that the `devtools` library is installed or install it by typing `install.packages("devtools")`
* Load the "devtools" library by typing 

`library(devtools)`

* Install the `HiddenSemiMarkov` package with the command:

`devtools::install_github("lasy/HiddenSemiMarkov")` or 

`devtools::install_github("lasy/HiddenSemiMarkov", build_vignettes = TRUE)` 
if you wish to build the vignette for this package.

* Find examples in the vignette by typing `browseVignettes(package = "HiddenSemiMarkov")` 
or start specifying a model with the `specify_hsmm` function (type `?specify_hsmm` for documentation).
