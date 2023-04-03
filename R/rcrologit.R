#' @title Estimation and Inference for Random Coefficients Rank-Ordered Logit
#'
#' @description This function takes the object prepared by \code{dataPrep} and estimates a 
#' random coefficient rank-ordered logit. The rank-ordered logit - sometimes termed \emph{exploded logit model} - was
#' originally proposed in \insertCite{beggs1981AssessingPotentialDemand;textual}{rcrologit} and it is an extension of the 
#' \insertCite{luce1959IndividualChoiceBehavior;textual}{rcrologit}-\insertCite{mcfadden1974FrontiersEconometrics;textual}{rcrologit} model. 
#'
#' These models automatically implies independence of irrelevant alternatives \insertCite{debreu1960ReviewRDLuce}{rcrologit}. However, 
#' \insertCite{mcfadden2000MixedMNLModels;textual}{rcrologit} show that if agents are allowed to have heterogeneous tastes (i.e., random coefficients), then the conditional
#' probability choices induced by the latent utility model can match those implied by virtually any
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
#'
#' @param dataPrep object of class 'rcrologit' prepared via \code{\link{dataPrep}}.
#' @param Sigma structure of the variance-covariance of the random coefficients. It must be one of
#' "diagonal" or "cholesky". Default is \code{Sigma="diagonal"}. See \strong{Details} section for more.
#' @param S integer denoting the number of simulations when approximating the integrals in the conditional
#' choice probabilities. Default is \code{S=50}.
#' @param approx.method string indicating the procedure to approximate the integrals in the conditional 
#' choice probabilities
#' @param stdErr string denoting whether standard error should be estimated and how. Available options are
#' "numerical" (default) which uses numerical approximation of the Hessian matrix; "analytical" which uses
#' the closed form of the variance of the score and the hessian of the likelihood to estimate standard
#' errors robust to misspecification; "skip" which makes R skip the computation of standard errors
#' @param verbose if \code{TRUE} prints additional information in the console. 
#' @param control.opts a list containing options to be passed to the underlying optimizer
#' \code{optim}.
#' 
#' @details 
#' 
#' \itemize{
#' \item{\strong{Variance-Covariance of random parameters.} The option \code{Sigma} allows the user to 
#' model directly the covariance structure of the random coefficients. Precisely, it shape \eqn{\boldsymbol{\Sigma}} in
#' \deqn{\left[\begin{array}{c}\boldsymbol{\alpha}_i \\ \boldsymbol{\beta}_i \end{array}\right]\sim
#' \mathsf{N}\left(\left[\begin{array}{c}\boldsymbol{\alpha}_{\mathtt{R}} \\ \boldsymbol{\beta}_{\mathtt{R}}\end{array}\right],
#' \boldsymbol{\Sigma}\right).}
#' In practice, to avoid issues with the positive-definiteness of \eqn{\boldsymbol{\Sigma}} during the optimization routine,
#' it is more convenient to model the random coefficients as
#' \deqn{\left[\begin{array}{c}\boldsymbol{\alpha}_i \\ \boldsymbol{\beta}_i \end{array}\right]=
#' \left[\begin{array}{c}\boldsymbol{\alpha}_{\mathtt{R}} \\ \boldsymbol{\beta}_{\mathtt{R}}\end{array}\right] + \Lambda u_i,
#' \quad u_i\sim\mathsf{N}(0,\boldsymbol{I}).}
#' The option \code{Sigma} directly models the shape of \eqn{\Lambda} and allows the user to choose between a diagonal structure
#' (\code{Sigma="diagonal"}) and a lower-triangular structure (\code{Sigma="cholesky"}).
#' }
#' }
#' 
#' @return
#' The function returns a list containing the following objects:
#' \item{b}{vector containing all the estimated parameters}
#' \item{bfix}{vector containing the estimated parameters corresponding to the "fixed" coefficients}
#' \item{bhet}{vector containing the estimated parameters corresponding to the "random" coefficients}
#' \item{Lambda}{loading matrix of the shocks to the random coefficients}
#' \item{Sigma}{variance-covariance matrix of the estimated parameters}
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

 
rcrologit <- function(dataPrep, Sigma="diagonal", S=50, approx.method="MC",
                       stdErr="numerical", verbose=FALSE, control.opts=NULL) {

  ################################################################################
  ## error checking
  if (methods::is(dataPrep, "dataRologit") == FALSE ) {
    stop("dataPrep should be the object returned by running dataPrep!")
  }
  
  if (!(Sigma %in% c("diagonal", "cholesky"))) {
    stop("The option 'Ltype' must be either 'diagonal' or 'cholesky'!")
  }
  
  if (!(approx.method %in% c("MC"))) stop("The option 'approx.method' must be either 'MC' or XXX")
  
  if (!(stdErr %in% c("skip", "analytical", "numerical"))) {
    stop("'stdErr' must be either 'skip', 'analytical', or 'numerical'!")
  }
  
  
  ################################################################################
  ## prepare list of matrices (X_1, X_2, ..., X_J)
  
  if (dataPrep$param.spec$K.het > 0) {
    rCoefs <- TRUE
    K.het.mu <- dataPrep$param.spec$K.het          # number of parameters following normal distribution (mean)
    if (Sigma=="diagonal") K.het.lam <- K.het.mu   # number of parameters in loading matrix 
    if (Sigma=="cholesky") K.het.lam <- K.het.mu*(K.het.mu+1)/2              
  } else {
    rCoefs <- FALSE
    K.het.mu <- K.het.lam <- 0
  }

  K.fix <- dataPrep$param.spec$K.fix  # number of fixed taste parameters to be estimated
  J <- dataPrep$param.spec$J          # total number of alternatives
  
  df <- data.frame("Worker.ID"=dataPrep$id,
                   "rank"=dataPrep$rank,
                   dataPrep$X.fix,
                   dataPrep$X.het)
  
  XX <- reshape(df, idvar = "Worker.ID",
                timevar = "rank",
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
  if (verbose) {
    cat("--------------------------------------------------------------------\n")
  }
  
  if (rCoefs == TRUE) { # random coefficients rank-ordered logit
    cat("Random Coefficients Rank-Ordered Logit \n")
    cat("Optimizing Likelihood function...\n")
    
    b0 <- rep(1, K.fix + K.het.mu + K.het.lam)     # order is going to be (beta.fix, beta.het, Sigma)
    
    if (approx.method == "MC") {     # draw shocks once and for all to approximate integrals
      epsMC <- MASS::mvrnorm(n=S, mu=rep(0, K.het.mu), Sigma=diag(K.het.mu))
    }
    
    #init <- Sys.time()
    bhat <- stats::optim(par=b0, fn=loglkldRC, X=Xlist, J=J, K.fix=K.fix, K.het.mu=K.het.mu, K.het.lam=K.het.lam,
                         Sigma=Sigma, approx.method="MC", S=S, epsMC=epsMC, method = "BFGS",
                         control=control.opts)
  
    ##########################################################
    # prepare output
    
    b <- bhat$par
    
    # unpack b to store results
    k1 <- K.fix                            
    k2 <- K.fix + 1
    k3 <- k1 + K.het.mu
    k4 <- k3 + 1
    k5 <- k3 + K.het.lam
    bfix <- b[1:k1]                                 # first sub-component is b for fixed taste covs (including FE)
    bhet <- b[k2:k3]                                # second: mean het taste 
    bLam <- vec2mat(b[k4:k5], K.het.mu, Sigma)      # third: varcov of shocks het taste 
    bLam <- abs(bLam)                               # loadings are identified up to a sign
   
    names(bfix) <- colnames(dataPrep$X.fix)
    names(bhet) <- colnames(dataPrep$X.het)
    rownames(bLam) <- colnames(bLam) <- paste0("Lambda.", names(bhet))  

  } else {  # rank-ordered logit
    cat("Rank-Ordered Logit \n")
    
    b0 <- rep(1, ncol(dataPrep$X.fix))
    bhat <- optim(par=b0, fn=loglkld, X=Xlist, method = "BFGS")
    b <- bhat$par
    names(b) <- colnames(dataPrep$X.fix)
    
    bfix <- b
    bhet <- bLam <- NULL
  }
  
  if (bhat$convergence!=0) {
    print(bhat$message)
    stop("Algorithm does not converge!")
  } else {
    cat("Algorithm reached convergence!\n")
  }
  
  if (verbose==TRUE) cat("Computing standard errors...\n\n")
  
  if (stdErr == "numerical") {  # numerical approximation
    
    se <- seGet(Xlist, b, rCoefs, NumApprox=TRUE, pars=NULL)
    SigmaHat <- se$Sigma

  } else if (stdErr == "analytical") { # analytical SEs
    se <- seGet(dataPrep, b, rCoefs, NumApprox=FALSE)
    SigmaHat <- se$SigmaRob
    
  } else {
    SigmaHat <- NULL

  }
  
  if (verbose==TRUE && stdErr!="skip") {
    bhat <- b
    sebhat <- sqrt(diag(SigmaHat))
    lb <- bhat - qnorm(0.975) * sebhat
    ub <- bhat + qnorm(0.975) * sebhat
    aux <- cbind(bhat, sebhat, lb, ub)
    colnames(aux) <- c("Estimate", "Std.Error", "Lb (95%)", "Ub (95%)")
    print(round(aux, digits = 2))
  }

  return(list(b=b, bfix=bfix, bhet=bhet, Lambda=bLam, Sigma=SigmaHat))
}
