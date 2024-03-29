# Overview
The purpose of PH2EWMA is to build a Phase II EWMA control chart with one-way random effects model in R.

# Installation

## Install from GitHub

PH2EWMA is still under development, so we recommend users installing it through Github as follows

``` r
install.packages("devtools")
devtools::install_github("bolus123/PH2EWMA")
```

## Install from local

Users can also download our release on Github and then install it from your local path as follows
``` r
install.packages('path_to_file', repos = NULL, type="source")
```


# Usage

PH2EWMA provides a function to build Phase II EWMA chart with variance components model as follows

``` r
data(Ph1data)
# Phase I sample

data(Ph2data)
# Phase II data

X1 <- as.matrix(Ph1data[, 2:4]) ^ (1/3)
X2 <- as.matrix(Ph2data[, 2:4]) ^ (1/3)
X2[which(is.na(X2))] <- mean(X1)

PH2EWMA(X2 = X2, X1 = X1) 
```

Also, PH2EWMA provides a function to get the corrected EWMA charting constant using the CUC or the EPC method as follows

``` r

ub <- c4.f(49)

# get the charting constant using the CUC method
getCC(m = 50, nu = 49, cc.option = 'CUC', ubCons = ub)

# get the charting constant using the EPC method
getCC(m = 50, nu = 49, cc.option = 'EPC', ubCons = ub)
```

More details are on the manual.
