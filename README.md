
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rcrologit

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-0.0.1-blue.svg)](https://github.com/filippopalomba/rcrologit)
[![R build
status](https://github.com/filippopalomba/rcrologit/workflows/R-CMD-check/badge.svg)](https://github.com/filippopalomba/rcrologit/actions)
[![](https://img.shields.io/github/last-commit/filippopalomba/rcrologit.svg)](https://github.com/filippopalomba/rcrologit/commits/master)
<!-- badges: end -->

The package provides estimation and inferential procedures for
rank-ordered logit model with agents with heterogeneous taste
preferences.

## Setup

The package estimates random coefficient rank-ordered logit models
induced by the latent utility model (McFadden, 1974) of the form
$$U_{i\ell} = u_{i\ell} + \epsilon_{i\ell},\quad i=1,2,\ldots,n,\quad j=0,1,\ldots,J.$$
In its most general form, we allow the user to model $u_{i\ell}$ as
$$u_{i\ell}=X_{i\ell}^\top\boldsymbol{\beta}_{\mathtt{F}} + Z_i^\top\boldsymbol{\alpha}_{\mathtt{F}} + 
 W_{i\ell}^\top\boldsymbol{\beta}_i + V_i^\top\boldsymbol{\alpha}_i + \delta_\ell$$
where - $X_{i\ell}$ are covariates varying at the unit-alternative level
whose coefficients are modelled as fixed - $Z_{i}$ are covariates
varying at the unit level whose coefficients are modelled as fixed -
$W_{i\ell}$ are covariates varying at the unit-alternative level whose
coefficients are modelled as random - $V_{i}$ are covariates varying at
the unit level whose coefficients are modelled as random - the random
coefficients are modeled as a joint multivariate normal and are i.i.d.
across units,
$$\left[\begin{array}{c}\boldsymbol{\alpha}_i \\ \boldsymbol{\beta}_i \end{array}\right]\sim
 \mathsf{N}\left(\left[\begin{array}{l}\boldsymbol{\alpha}_{\mathtt{R}} \\ \boldsymbol{\beta}_{\mathtt{R}}\end{array}\right],\boldsymbol{\Sigma}\right) $$ -
$\delta_\ell$ are alternative-specific fixed effects -
$\epsilon_{i\ell}\sim\mathrm{Gu}(0,1)}$ are idiosyncratic i.i.d. shocks

Note that whenever $W_{i\ell}$ and $V_i$ are not specified estimates a
standard rank-ordered logit with no heterogeneous preferences.

## Installation

You can install the development version of rcrologit from
[GitHub](https://github.com/filippopalomba/rcrologit) with:

``` r
# install.packages("devtools")
devtools::install_github("filippopalomba/rcrologit")
```

## Basic Usage

``` r
library(rcrologit)

data <- rcrologit_data

# Rank-ordered logit
dataprep <- dataPrep(data, idVar = "Worker_ID", rankVar = "rank",
                    altVar = "alternative",
                    covsInt.fix = list("Gender"),
                    covs.fix = list("log_Wage"), FE = c("Firm_ID"))
                      
rologitEst <- rcrologit(dataprep)

# Rank-ordered logit
dataprep <- dataPrep(data, idVar = "Worker_ID", rankVar = "rank",
                    altVar = "alternative",
                    covsInt.het = list("Gender"),
                    covs.fix = list("log_Wage"), FE = c("Firm_ID"))
                      
rologitEst <- rcrologit(dataprep, stdErr="skip")
```
