#' Simulated Data According to McFadden (1974) Latent Utility Model
#'
#' A dataset containing three-option rankings of 3000 units and some covariates
#'
#' @docType data
#'
#' @usage data(rcrologit_data)
#'
#' @format A data frame with 9000 rows and 11 variables:
#' \describe{
#'   \item{Worker_ID}{worker id.}
#'   \item{alternative}{job offer id (within worker).}
#'   \item{rank}{rank of the job offer.}
#'   \item{Firm_ID}{id of the firm offering the job.}
#'   \item{Gender}{gender of worker (1 if Female).}
#'   \item{Educ}{education level of worker (1 if <HS, 2 if HS, 3 if College or more).}
#'   \item{WorkerSector}{sector of worker's previous occupation (1 if Manufacturing, 0 otherwise).}
#'   \item{Wage}{wage offered.}
#'   \item{log_Wage}{log of wage offered.}
#'   \item{FirmSector}{sector of the firm offering the job (1 if Manufacturing, 0 otherwise).}
#' }
#' 
#' @keywords datasets
#' 					
"rcrologit_data"