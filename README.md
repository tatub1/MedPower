
# MedPower

<!-- badges: start -->
<!-- badges: end -->

The goal of MedPower is to calculate power or sample size to detect natural indirect effects for causal mediation analysis.

## Installation

You can install the development version of MedPower like so:

``` r
library(devtools)
devtools::install_github("tatub1/MedPower")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MedPower)
## basic example code

MedPower_CC(
  beta1 = sqrt(0.25),
  theta2 = 0.29,
  theta3 = 0,
  power = 0.8
)
```

