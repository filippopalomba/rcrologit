#' @title \code{rcrologit}: A Package to Estimate Random Coefficients Rank-Ordered Logit Models.
#'
#' @description The package implements estimation and inference procedures for random coefficient rank-ordered logit models.
#' 
#' @author 
#' Chiara Motta, University of California Berkeley. \email{cmotta@berkeley.edu}
#' 
#' Filippo Palomba, Princeton University. \email{fpalomba@princeton.edu}
#'
#' @importFrom fastDummies dummy_cols
#' @importFrom MASS mvrnorm
#' @importFrom methods is
#' @importFrom numDeriv hessian
#' @importFrom parallel mclapply
#' @importFrom purrr map
#' @importFrom Rdpack reprompt
#' @importFrom tibble is_tibble
#' 
#' @rawNamespace import(stats, except = c(lag, filter, power))
#' @rawNamespace import(rlang, except = c(is_vector, is_complex))
#'
#' @docType package
#'
#' @aliases rcrologit-package
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1")  utils::globalVariables(c("."))