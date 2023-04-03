#' @title Data Preparation for Estimation of Random Coefficients Rank-Ordered Logit
#'
#' @description This function prepares the data for estimation of a random coefficient rank-ordered
#' logit. The rank-ordered logit - sometimes termed \emph{exploded logit model} - was
#' originally proposed in \insertCite{beggs1981AssessingPotentialDemand;textual}{rcrologit} and it is an extension of the
#' \insertCite{luce1959IndividualChoiceBehavior;textual}{rcrologit}-\insertCite{mcfadden1974FrontiersEconometrics;textual}{rcrologit} model. 
#'
#' These models automatically implies independence of irrelevant alternatives \insertCite{debreu1960ReviewRDLuce}{rcrologit}. However, 
#' \insertCite{mcfadden2000MixedMNLModels;textual}{rcrologit} show that if agents are allowed to have heterogeneous tastes
#' (i.e., random coefficients),
#' then the conditional probability choices induced by the latent utility model can match those implied by virtually any
#' discrete choice probability model.
#'
#' The package \code{\link{rcrologit}}, depending on the type of covariates specified in \code{dataPrep},
#' allows the user to estimate:
#' \itemize{
#' \item{the standard rank-ordered logit model if either \code{covs.fix} or \code{covsInt.fix} are specified.}
#' \item{the random coefficients rank-ordered logit model if either \code{covs.het} or \code{covsInt.het} are specified.}
#' \item{the mixed random coefficients rank-ordered logit model if at least one of \code{covs.fix} or \code{covsInt.fix} and
#' at least one of \code{covs.het} or \code{covsInt.het} are specified.}
#' }
#'
#' For more information on the underlying specification see the \href{https://github.com/filippopalomba/rcrologit}{official repository}.
#'
#' @param dataraw a \code{data.frame} object containing the data to be prepared for estimation
#' @param idVar a string indicating the name of the unit identifier variable
#' @param rankVar a string indicating the name of the variable indicating the rank of each alternative.
#' It should be a variable containing integer values ranging from 1 to \eqn{J}, where \eqn{J} is the number
#' of alternatives. Each unit's best alternative must get value 1, the second best value 2, and so on in a decreasing
#' fashion.
#' @param altVar a string indicating the name of the alternative identifier variable
#' @param FE a string containing the name of the fixed effect variable (if required by the user)
#' @param covs.fix covariates varying at the alternative-unit level with fixed taste parameter
#' @param covs.het covariates varying at the alternative-unit level with random taste parameter
#' @param covsInt.fix covariates varying at the unit level with fixed taste parameter
#' @param covsInt.het covariates varying at the unit level with random taste parameter
#'
#' @details
#' The function prepares the data to estimate a random coefficient rank-ordered logit model induced by the latent utility model 
#' \insertCite{mcfadden1974FrontiersEconometrics;textual}{rcrologit}
#' \deqn{U_{i\ell} = u_{i\ell} + \epsilon_{i\ell},\quad i=1,2,\ldots,n,\quad j=0,1,\ldots,J.}
#' In its most general form, we model \eqn{u_{i\ell}} as
#' \deqn{u_{i\ell}=X_{i\ell}^\top\boldsymbol{\beta}_{\mathtt{F}} + Z_i^\top\boldsymbol{\alpha}_{\ell,\mathtt{F}} +
#' W_{i\ell}^\top\boldsymbol{\beta}_i + V_i^\top\boldsymbol{\alpha}_{i\ell} + \delta_\ell}
#' where
#' \itemize{
#' \item{\eqn{X_{i\ell}} are covariates varying at the unit-alternative level whose coefficients are modelled as fixed.
#' The user can specify these covariates via the option \code{covs.fix}.}
#' \item{\eqn{Z_{i}} are covariates varying at the unit level whose coefficients are modelled as fixed.
#' The user can specify these covariates via the option \code{covsInt.fix}. These coefficients are interacted with
#' a dummy for the choice and treated as alternative varying covariates,
#' i.e., \eqn{Z_{i\ell}=\sum_{j=1}^JZ_i\times\mathbf{1}(j=\ell)}, where $J=0$ is the reference group.}
#' \item{\eqn{W_{i\ell}} are covariates varying at the unit-alternative level whose coefficients are modelled as random.
#' The user can specify these covariates via the option \code{covs.het}.}
#' \item{\eqn{V_{i}} are covariates varying at the unit level whose coefficients are modelled as random.
#' The user can specify these covariates via the option \code{covsInt.het}. These coefficients are interacted with
#' a dummy for the choice and treated as alternative varying covariates,
#' i.e., \eqn{V_{i\ell}=\sum_{j=1}^JV_i\times\mathbf{1}(j=\ell)}, where $J=0$ is the reference group.}
#' \item{the random coefficients}{ are modeled as a joint multivariate normal and are i.i.d. across units,
#' \deqn{\left[\begin{array}{c}\boldsymbol{\alpha}_i \\ \boldsymbol{\beta}_i \end{array}\right]\sim
#' \mathsf{N}\left(\left[\begin{array}{l}\boldsymbol{\alpha}_{\mathtt{R}} \\
#' \boldsymbol{\beta}_{\mathtt{R}}\end{array}\right],\boldsymbol{\Sigma}\right) }}
#' \item{\eqn{\delta_\ell}}{ are alternative-specific fixed effects that can be specified via the option \code{FE}.}
#' \item{\eqn{\epsilon_{i\ell}\sim\mathsf{Gu}(0,1)}}{are idiosyncratic i.i.d. shocks.}
#' }
#' The parameter vector to be estimated is thus
#' \deqn{\theta = \left(\boldsymbol{\beta_\mathtt{F}}^\top,\boldsymbol{\beta_\mathtt{R}}^\top,
#' \boldsymbol{\alpha_\mathtt{F}}^\top,\boldsymbol{\alpha_\mathtt{R}}^\top,
#' \mathrm{vech}(\boldsymbol{\Sigma})^\top,\{\delta\}_{j=1}^J\right)^\top,}
#' where the first alternative-fixed effect has been normalized to 0.
#'
#' For more information on the underlying specification see the \href{https://github.com/filippopalomba/rcrologit}{official repository}.
#'
#' @return
#' The function returns a list containing the following objects:
#' \item{id}{vector containing the id of each observation}
#' \item{rank}{vector containing the rank of each observation}
#' \item{X.fix}{a matrix containing covariates with fixed taste parameter}
#' \item{X.het}{a matrix containing covariates with random taste parameter}
#' \item{param.spec}{a list containing some parameters describing the specification chosen by the user}
#'
#' @author
#' Chiara Motta, University of California Berkeley. \email{cmotta@berkeley.edu}
#'
#' Filippo Palomba, Princeton University. \email{fpalomba@princeton.edu}
#'
#' @references
#'  \insertAllCited{}
#'
#' @seealso \code{\link{rcrologit}}
#'
#' @examples
#' data <- rcrologit_data
#' dataprep <- dataPrep(data, idVar = "Worker_ID", rankVar = "rank",
#'                      altVar = "alternative",
#'                      covsInt.fix = list("Gender"),
#'                      covs.fix = list("log_Wage"), FE = c("Firm_ID"))
#'
#' rologitEst <- rcrologit(dataprep)
#'
#' @export
#'
dataPrep <- function(dataraw, idVar, rankVar, altVar, FE=NULL,
                     covs.fix=NULL, covs.het=NULL, covsInt.fix=NULL, covsInt.het=NULL) {

  # copy dataset to avoid overwriting
  data <- dataraw

  # sanity checks that options are well-specified
  errorCheck(data, idVar, rankVar, altVar, FE, covs.fix, covs.het, covsInt.fix, covsInt.het)

  # make sure that data is sorted by unit and by rank (1,2,3,...)
  data <- data[order(data[[idVar]], data[[rankVar]]), ]

  # data preparation
  # select id and rank variable based on input names
  id <- data[, idVar]            # unit id
  rank <- data[, rankVar]        # rank

  # select covariates based on input name
  if (!is.null(covs.fix)) {X.fix <- as.matrix(data[, unlist(covs.fix)])} else {X.fix <- NULL}
  if (!is.null(covs.het)) {X.het <- as.matrix(data[, unlist(covs.het)])} else {X.het <- NULL}
  if (!is.null(covsInt.fix)) {X.fix.int <- as.matrix(data[, unlist(covsInt.fix)])} else {X.fix.int <- NULL}
  if (!is.null(covsInt.het)) {X.het.int <- as.matrix(data[, unlist(covsInt.het)])} else {X.het.int <- NULL}

  # add interaction between dummies for alternatives with covsInt (id)
  if (!is.null(covsInt.fix)) {
    D <- fastDummies::dummy_cols(data,
                                 select_columns = altVar,
                                 remove_first_dummy = TRUE)
    D <- as.matrix(D[,(ncol(data)+1):ncol(D)])
    colnames(D) <- c(2:max(rank))

    Xint <- c()
    for (var in covsInt.fix) {
      aux <- D*data[[var]]
      colnames(aux) <- paste(var, colnames(D), sep="X")
      Xint <- cbind(Xint, aux)
    }
    X.fix <- cbind(X.fix, Xint)
  }

  if (!is.null(covsInt.het)) {
    D <- fastDummies::dummy_cols(data,
                                 select_columns = altVar,
                                 remove_first_dummy = TRUE)
    D <- as.matrix(D[,(ncol(data)+1):ncol(D)])
    colnames(D) <- c(2:max(rank))

    Xint <- c()
    for (var in covsInt.het) {
      aux <- D*data[[var]]
      colnames(aux) <- paste(var, colnames(D), sep="X")
      Xint <- cbind(Xint, aux)
    }
    X.het <- cbind(X.het, Xint)
  }

  # create fixed effects
  if (!is.null(FE)) {
    dataFE <- fastDummies::dummy_cols(data,
                                      select_columns = FE,
                                      remove_first_dummy = TRUE)
    Xfe <- dataFE[,c((ncol(data)+1): ncol(dataFE))]
    
    if (!is.null(X.fix)) {X.fix <- cbind(X.fix, Xfe)} else {X.fix <- Xfe}
    Kfe <- ncol(Xfe) 
  } else {
    Kfe <- 0
  }
  
  J <- length(unique(data[, altVar])) # number of alternatives
  K_il.fix <- length(covs.fix)        # number of unit-alternative varying covs with fixed taste
  K_il.het <- length(covs.het)        # number of unit-alternative varying covs with het taste
  K_i.fix <- length(covsInt.fix)      # number of unit varying covs with fixed taste
  K_i.het <- length(covsInt.het)      # number of unit varying covs with het taste
  K <- (K_il.fix + K_il.het) + (J-1)*(K_i.fix + K_i.het) # total number of covs
  
  # store a bunch of useful parameters
  param.spec <-   list(K_il.fix = K_il.fix,                    # number of unit-alternative varying covs with fixed taste
                       K_il.het = K_il.het,                    # number of unit-alternative varying covs with het taste
                       K_i.fix = K_i.fix,                      # number of unit varying covs with fixed taste
                       K_i.het = K_i.het,                      # number of unit varying covs with het taste
                       K.fix = K_il.fix + (J-1)*K_i.fix + Kfe, # total number of covs with fixed taste
                       K.het = K_il.het + (J-1)*K_i.het,       # total number of covs with het taste
                       K = K,                                  # total number of covs
                       J = J,                                  # number of alternatives 
                       N = length(unique(data[, idVar])))        # number of observations

  if (!is.null(X.fix)) {
    X.fix <- as.matrix(X.fix)
  } else {
    X.fix <- data.frame(matrix(nrow = param.spec$N*J, ncol = 0)) 
  }
  if (!is.null(X.het)) {
    X.het <- as.matrix(X.het)
  } else {
    X.het <- data.frame(matrix(nrow = param.spec$N*J, ncol = 0))
  }
  toret <- list(id = id, rank = rank, X.fix = X.fix,
                X.het = X.het, param.spec=param.spec)
  
  class(toret) <- 'dataRologit'
  return(toret)
}
