
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

This is a basic example which shows you how to solve a sample-size problem:

``` r
library(MedPower)
## Basic sample-size calculation example code
## We assume mediator and outcome are continuous.
## default values are set as all variances=1, alpha=0.05.
## You have to specify either power or sample size.

MedPower_CC(
  beta1 = sqrt(0.25),
  theta2 = 0.29,
  theta3 = 0,
  power = 0.8
)
```

