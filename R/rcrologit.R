#' @title Estimation and Inference for Random Coefficients Rank-Ordered Logit
#'
#' @description This function takes the object prepared by \code{dataPrep} and estimates a
#' random coefficient rank-ordered logit. The rank-ordered logit - sometimes termed \emph{exploded logit model} - was
#' originally proposed in \insertCite{beggs1981AssessingPotentialDemand;textual}{rcrologit} and it is an extension of the
#' \insertCite{luce1959IndividualChoiceBehavior;textual}{rcrologit}-\insertCite{mcfadden1974FrontiersEconometrics;textual}{rcrologit} model. 
#'
#' These models automatically implies independence of irrelevant alternatives \insertCite{debreu1960ReviewRDLuce}{rcrologit}. However, 
#' \insertCite{mcfadden2000MixedMNLModels;textual}{rcrologit} show that if agents are allowed to have heterogeneous
#' tastes (i.e., random coefficients), then the conditional probability choices induced by the latent utility model can
#' match those implied by virtually any discrete choice probability model.
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
#'
#' @param dataprep object of class 'rcrologit' prepared via \code{\link{dataPrep}}.
#' @param Sigma structure of the variance-covariance of the random coefficients. It must be one of
#' "diagonal" or "cholesky". Default is \code{Sigma="diagonal"}. See \strong{Details} section for more.
#' @param S integer denoting the number of simulations when approximating the integrals in the conditional
#' choice probabilities. Default is \code{S=50}. More on how to select the optimal \eqn{S} can be found in
#' \insertCite{cameron2005MicroeconometricsMethodsApplications}{rcrologit} and \insertCite{hajivassiliou2000PracticalIssuesMaximum}{rcrologit}.
#' @param approx.method string indicating the procedure to approximate the integrals in the conditional
#' choice probabilities when including random coefficients.
#' At the moment only approximation via monte-carlo simulation is available. Future releases will include
#' importance sampling and other alternatives
#' @param bias.correction if \code{TRUE} applies the bias correction for simulated maximum likelihood proposed in 
#' \insertCite{gourieroux1991SimulationBasedInference}{rcrologit}. For more details see Section 12.4.4 in \insertCite{cameron2005MicroeconometricsMethodsApplications}{rcrologit}.
#' This option is effective only when random coefficients are 
#' included in the model.
#' @param robust if \code{TRUE} computes standard errors robust to misspecification.
#' @param stdErr.dfadj boolean indicating whether a degrees-of-freedom correction should be used when estimating
#' the variance of the ML estimator. See \strong{Details} section for more. 
#' @param skip.stdErr if \code{TRUE} skips computation of standard errors which can be intensive in the random 
#' coefficient rologit model.
#' @param Ncores integer indicating the number of cores to be used in simulating conditional choice probabilities.
#' It affects speed only when random coefficients are included in the model. Speed gains are sensible whenever \eqn{S\geq 100}.
#' On Windows systems it is set automatically to \code{Ncores = 1}.
#' @param verbose if \code{TRUE} prints additional information in the console.
#' @param control.opts a list containing options to be passed to the underlying optimizer \code{optim}.
#'
#' @details
#' \itemize{
#' \item{\strong{Variance-Covariance degrees-of-freedom adjustment.} When \code{stdErr.dfadj = FALSE}, the estimate for the 
#' asymptotic variance of \eqn{\sqrt{N}(\widehat{\boldsymbol{\theta}}-\boldsymbol{\theta})} is divided by the sample size \eqn{N}. If instead,
#' \code{stdErr.dfadj = TRUE}, then the estimate for the 
#' asymptotic variance of \eqn{\sqrt{N}(\widehat{\boldsymbol{\theta}}-\boldsymbol{\theta})} is divided by \eqn{N-k}, where \eqn{k} is the 
#' dimension of \eqn{\boldsymbol{\theta}}.}
#' 
#' \item{\strong{Variance-Covariance of random parameters.} The option \code{Sigma} allows the user to
#' model directly the covariance structure of the random coefficients. Precisely, it shape \eqn{\boldsymbol{\Sigma}} in
#' \deqn{\left[\begin{array}{c}\boldsymbol{\alpha}_i \\ \boldsymbol{\beta}_i \end{array}\right]\sim
#' \mathsf{N}\left(\left[\begin{array}{c}\boldsymbol{\alpha}_{\mathtt{R}} \\ \boldsymbol{\beta}_{\mathtt{R}}\end{array}\right],
#' \boldsymbol{\Sigma}\right).}
#' In practice, to avoid issues with the positive-definiteness of \eqn{\boldsymbol{\Sigma}} during the optimization routine,
#' it is more convenient to model the random coefficients as
#' \deqn{\left[\begin{array}{c}\boldsymbol{\alpha}_i \\ \boldsymbol{\beta}_i \end{array}\right]=
#' \left[\begin{array}{c}\boldsymbol{\alpha}_{\mathtt{R}} \\ \boldsymbol{\beta}_{\mathtt{R}}\end{array}\right] + \Lambda \mathbf{u}_i,
#' \quad \mathbf{u}_i\sim\mathsf{N}(0,\boldsymbol{I}).}
#' The option \code{Sigma} directly models the shape of \eqn{\Lambda} and allows the user to choose between a diagonal structure
#' (\code{Sigma="diagonal"}) and a lower-triangular structure (\code{Sigma="cholesky"}).
#' }
#' \item{\strong{Model.}
#' The function estimates a random coefficient rank-ordered logit model induced by the latent utility model
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
#' where the first alternative-fixed effect has been normalized to 0.}
#' }
#'
#' For more information on the underlying specification see the \href{https://github.com/filippopalomba/rcrologit}{official repository}.
#'
#' @return
#' The function returns a list containing the following objects:
#' \item{b}{vector containing all the estimated parameters}
#' \item{bfix}{vector containing the estimated parameters corresponding to the "fixed" coefficients}
#' \item{bhet}{vector containing the estimated parameters corresponding to the "random" coefficients}
#' \item{Lambda}{loading matrix of the shocks to the random coefficients}
#' \item{Sigma}{variance-covariance matrix of the estimated parameters}
#' \item{param.spec}{a list containing some parameters describing the specification chosen by the user}
#'
#' @author
#' Chiara Motta, University of California Berkeley. \email{cmotta@berkeley.edu}
#'
#' Filippo Palomba, Princeton University. \email{fpalomba@princeton.edu}
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
#' @seealso \code{\link{dataPrep}}
#'
#' @references
#'  \insertAllCited{}
#'
#' @export

rcrologit <- function(dataprep, Sigma = "diagonal", S = 50L, approx.method = "MC", bias.correction = FALSE,
                      robust = FALSE, stdErr.dfadj = TRUE, skip.stdErr = FALSE, Ncores = 1L,
                      verbose = TRUE, control.opts = NULL) {

  ################################################################################
  ## error checking
  if (methods::is(dataprep, "dataRologit") == FALSE) {
    stop("dataprep should be the object returned by running dataprep!")
  }

  if (!(Sigma %in% c("diagonal", "cholesky"))) {
    stop("The option 'Ltype' must be either 'diagonal' or 'cholesky'!")
  }

  if (!(approx.method %in% c("MC"))) stop("The option 'approx.method' must be either 'MC' or XXX")

  if (.Platform$OS.type != "unix") Ncores <- 1
  
  ################################################################################
  ## prepare list of matrices (X_1, X_2, ..., X_J)
  rCoefs <- dataprep$param.spec$model == "rcoef_rologit"
  if (rCoefs == TRUE) {
    K.het.mu <- dataprep$param.spec$K.het          # number of parameters following normal distribution (mean)
    if (Sigma == "diagonal") K.het.lam <- K.het.mu   # number of parameters in loading matrix
    if (Sigma == "cholesky") K.het.lam <- K.het.mu * (K.het.mu + 1) / 2
  } else {
    K.het.mu <- K.het.lam <- 0
  }

  K.fix <- dataprep$param.spec$K.fix  # number of fixed taste parameters to be estimated
  J <- dataprep$param.spec$J          # total number of alternatives

  df <- data.frame(dataprep$id,
                   dataprep$rank,
                   dataprep$X.fix,
                   dataprep$X.het)
  colnames(df) <- c(dataprep$param.spec$idVar, dataprep$param.spec$rankVar, colnames(df[3:ncol(df)]))

  XX <- reshape(df, idvar = dataprep$param.spec$idVar,
                timevar = dataprep$param.spec$rankVar,
                direction = "wide")
  
  XX <- as.matrix(XX[,-1])
  m <- K.fix + K.het.mu 
  
  Xlist <- list()
  for (j in seq_len(J)) {
    jl <- (j-1)*m + 1
    ju <- j*m
    Xlist[[j]] <- XX[,jl:ju]
  }  
  
  ################################################################################
  ## Maximum Likelihood Estimation
  if (verbose) cat("--------------------------------------------------------------------\n")
  

  ########################################
  # random coefficients rank-ordered logit
  ########################################

  if (rCoefs == TRUE) { 
    if (verbose) cat("Random Coefficients Rank-Ordered Logit \n")
    if (verbose) cat("Optimizing Likelihood function...\n")

    b0 <- rep(1, K.fix + K.het.mu + K.het.lam)     # order is going to be (beta.fix, beta.het, Sigma)

    if (approx.method == "MC") {     # draw shocks once and for all to approximate integrals
      epsMC <- MASS::mvrnorm(n=S, mu=rep(0, K.het.mu), Sigma=diag(K.het.mu))
    }

    bhat <- stats::optim(par=b0, fn=loglkldRC, X=Xlist, J=J, K.fix=K.fix, K.het.mu=K.het.mu, K.het.lam=K.het.lam,
                         Sigma=Sigma, bias.corr=bias.correction, approx.method="MC", S=S, epsMC=epsMC, Ncores=Ncores,
                         method = "BFGS", control=control.opts)
    
    # prepare output
    b <- as.matrix(bhat$par)

    # unpack b to store results
    k1 <- K.fix
    k2 <- K.fix + 1
    k3 <- k1 + K.het.mu
    k4 <- k3 + 1
    k5 <- k3 + K.het.lam
    bfix <- b[1:k1, drop=FALSE]                                 # first sub-component is b for fixed taste covs (including FE)
    bhet <- b[k2:k3, drop=FALSE]                                # second: mean het taste
    bLam <- vech2mat(b[k4:k5, drop=FALSE], K.het.mu, Sigma)     # third: varcov of shocks het taste
    bLam <- as.matrix(abs(bLam))                                # loadings are identified up to a sign

    names(bfix) <- colnames(dataprep$X.fix)
    names(bhet) <- colnames(dataprep$X.het)
    rownames(bLam) <- colnames(bLam) <- paste0("Lambda.", names(bhet))
  }
  
  ########################################
  # rank-ordered logit
  ########################################
    
  if (rCoefs == FALSE) {  
    if (verbose) cat("Rank-Ordered Logit \n")
    
    b0 <- rep(1, ncol(dataprep$X.fix))
    bhat <- stats::optim(par=b0, fn=loglkld, X=Xlist, J=J, method = "BFGS")
    b <- bhat$par
    names(b) <- colnames(dataprep$X.fix)
    
    bfix <- b
    bhet <- bLam <- NULL
    
    # retrieve fitted values
    #browser()
  }

  
  ################################################################################
  ## Convergence check
  
  if (bhat$convergence!=0) {
    print(bhat$message)
    stop("Algorithm has not converged!")
  } else {
    if (verbose) cat("Algorithm reached convergence!\n")
  }
  
  ################################################################################
  ## Standard Error Computation
  
  if (skip.stdErr == TRUE) {
    
    SigmaHat <- matrix(NA, nrow(b), nrow(b))
    
  } else { 
    
    if (verbose) cat("Computing standard errors...\n")
    
    ########################################
    # rank-ordered logit
    ########################################    
    
    if (rCoefs == FALSE) {
      
      pars <- list(N = dataprep$param.spec$N, J = dataprep$param.spec$J)
      se <- seGet(Xlist, b, rCoefs, robust=robust, pars=pars)
      SigmaHat <- se$Sigma
      
    }
    
    ########################################
    # random coefficients rank-ordered logit
    ########################################    
    
    if (rCoefs == TRUE) {
      
      pars <- list(bfix=bfix, bhet=bhet, bLam=bLam, J=J, bias.corr=bias.correction, 
                   approx.method=approx.method, S=S, epsMC=epsMC, Ncores=Ncores,
                   K.fix=K.fix, K.het.mu=K.het.mu, K.het.lam=K.het.lam, Sigma=Sigma,
                   N=dataprep$param.spec$N, verbose=verbose)
      se <- seGet(Xlist, b, rCoefs, robust=robust, pars=pars)
      SigmaHat <- se$Sigma
      
    }
    
    if (stdErr.dfadj == TRUE ) {
      N <- dataprep$param.spec$N
      k <- length(b)
      SigmaHat <- SigmaHat * ((N - k) / N)
    }
  }

  ################################################################################
  ## Return stuff
  
  to_return <- list(b = b, bfix = bfix, bhet = bhet,
                    Lambda = bLam, Sigma = SigmaHat, param.spec = dataprep$param.spec)
  to_return$param.spec$K.het.lam <- K.het.lam
  to_return$param.spec$Sigma <- Sigma
  to_return$param.spec$robust <- robust

  class(to_return) <- "rcrologit"

  return(to_return)
}
