##########################################################################
##########################################################################
# computes the likelihood of the rank-ordered logit
loglkld <- function(b0, X) {
  
  eXb <- lapply(X, function(x) exp(x%*%b0))
  lkld <- log(eXb[[1]]/(eXb[[1]] + eXb[[2]] + eXb[[3]])) + log(eXb[[2]]/(eXb[[2]] + eXb[[3]]))
  
  return(-sum(lkld))
}

##########################################################################
# computes the likelihood of the random coefficients rank-ordered logit
loglkldRC <- function(b0, X, J, K.fix, K.het.mu, K.het.lam, Sigma, approx.method, S, epsMC) {

  # unpack b0 to pass it to mvnorm
  k1 <- K.fix                            
  k2 <- K.fix + 1
  k3 <- k1 + K.het.mu
  k4 <- k3 + 1
  k5 <- k3 + K.het.lam
  bfix <- b0[1:k1]                                # first sub-component is b for fixed taste covs (including FE)
  bhet <- b0[k2:k3]                               # second: mean het taste 
  bLam <- vech2mat(b0[k4:k5], K.het.mu, Sigma)    # third: varcov of shocks het taste 

  ccp <- ccpGet(bfix, bhet, bLam, X, J, approx.method, S, epsMC) 
  
  lkld <- -sum(log(ccp))   
  
  return(lkld)
}

##########################################################################
# estimates the conditional choice probabilities
ccpGet <- function(bfix, bhet, bLam, X, J, approx.method, S, epsMC) {
  
  if (approx.method == "MC") {
    
    # shift and rescale normal iid draws
    b_i <- t(apply(epsMC, 1, function(x) bLam %*% x))  # apply loading matrix to each row of shocks
    b_i <- sweep(b_i, 2, bhet, "+")                    # add new mean to scaled shocks
    if (nrow(b_i) == 1) b_i <- t(b_i)
    
    # replicate bfix over each row to make it compatible with b_i
    bfix <- matrix(bfix, nrow=S, ncol=length(bfix), byrow=TRUE)
    b_i <- cbind(bfix, b_i)  # fixed taste parameters come first
    
    # apply to each row (i.e. draw) the ccpROLogit function to compute ccp of
    # each unit's ranking. apply is going to store the ccp in a column vector
    # and then cbinds all of then. thus to get the MC approximation of the 
    # integral we average the result over columns (i.e. draws)
    ccp <- rowMeans(apply(b_i, 1, function(x) ccpROLogit(X, x, J)))
  }
  
  return(ccp)
}

##########################################################################
# compute rank ordered logit ccp
ccpROLogit <- function(X, b, J) {
  eXb <- lapply(X, function(x) exp(x%*%as.matrix(b)))
  rol <- (eXb[[1]]/(eXb[[1]] + eXb[[2]] + eXb[[3]])) * (eXb[[2]]/(eXb[[2]] + eXb[[3]]))
  return(rol)
}

##########################################################################
# since mat will be symmetric we don't care about the ordering (by row or by col)
vech2mat <- function(vec, dmn, shape) { 
  if (shape == "diagonal") {
    if (dmn == 1) {
      mat <- diag(as.matrix(vec))
    } else {
      mat <- diag(vec)
    }
  } else {
    m1 <- matrix(NA, dmn, dmn)
    m1[lower.tri(m1, diag=TRUE)] <- vec
    mat <- pmax(m1, t(m1), na.rm=TRUE)
  }
  return(mat)
}



##########################################################################
##########################################################################
# analytic standard errors
seGet <- function(Xlist, b, rCoefs, pars, NumApprox=TRUE) {
  
  if (NumApprox == TRUE) { # numerical approximation of hessian of likelihood function
    
    if (rCoefs == FALSE) {
      
      Hes <- numDeriv::hessian(loglkld, x = b,
                               X = Xlist)    
      Sigma <- solve(Hes)
      
    } else {
      
     id <- pars$id 
     rank <- pars$rank
     K <- pars$K
     Klam <- pars$Klam 
     Kfe <- pars$Kfe 
     approx.method <- pars$approx.method 
     S <- pars$S
      
      loglkld2 <- function(b0, X) {
        
        # unpack b0 to pass it to mvnorm
        k1 <- K
        k2 <- K+1
        k3 <- k1 + Klam
        k4 <- k3 + 1
        k5 <- k3 + Kfe
        bmean <- b0[1:k1]
        bLam <- vech2mat(b0[k2:k3], K)
        bfe <- b0[k4:k5]
        
        ccp <- ccpGet(bmean, bLam, bfe, X, id, rank, approx.method, S) 
        
        lkld <- -sum(log(ccp))   
        return(lkld)
      }
    }
    
    return(list(Sigma = Sigma, SigmaRob = NULL, G = NULL, Hes = Hes))
    
  } else {  # analytical formula
    
    id <- dataPrep$id
    idlist <- unique(id[[1]])
    rank <- dataPrep$rank 
    X <- dataPrep$X
    
    Xb <- X %*% b
    eXb <- exp(Xb)
    
    XX <- data.frame(id=id, rank=rank, X=X, eXb = eXb)
    colnames(XX) <- c("id", "rank", colnames(XX[,3:ncol(XX)]))
    
    Jac <- matrix(NA, length(idlist), ncol(X))
    Hes <- matrix(0, ncol(X), ncol(X))
    
    j <- 1
    for (i in idlist) {
      dataID <- subset(XX, id == i)
      
      # Computing variance of score V (meat)
      exb <- dataID$eXb
      eXbi <- sum(exb)
      eXbi2 <- sum(exb[2:3])
      Xi <- t(dataID[,3:(ncol(dataID)-1)])
      
      eXbX <- as.matrix(rowSums(Xi * exb))
      eXbX2 <- as.matrix(rowSums(Xi[,2:3] * exb[2:3]))
      
      aux <- (Xi[, 1] - eXbX / eXbi) + (Xi[, 2] - eXbX2 / eXbi2)
      Jac[j, ] <- t(aux)
      
      # Computing derivative of the score/Hessian (bread)
      h1 <- exb[1] * Xi[,1] %*% t(Xi[,1]) + exb[2] * Xi[,2] %*% t(Xi[,2]) + exb[3] * Xi[,3] %*% t(Xi[,3])
      h2 <- exb[2] * Xi[,2] %*% t(Xi[,2]) + exb[3] * Xi[,3] %*% t(Xi[,3])
      
      aux <- (eXbX %*% t(eXbX) - eXbi * h1 ) / eXbi^2 +  (eXbX2 %*% t(eXbX2) - eXbi2 * h2 ) / eXbi2^2
      Hes <- Hes + aux
      
      j <- j + 1
    }
    
    Hinv <- solve(Hes) # take average and get inverse of hessian
    G <- t(Jac) %*% Jac 
   
    return(list(Sigma = G / length(idlist), SigmaRob = Hinv %*% G %*% Hinv, G = G, Hes = Hes)) 
  }
}


##########################################################################
##########################################################################
# error checking in dataPrep
errorCheck <- function(data, idVar, rankVar, altVar, FE, covs.fix, covs.het, covsInt.fix, covsInt.het) {
  if (!is.character(idVar)) stop("idVar should be a character! (eg. idVar = 'ID')")
  if (!is.character(rankVar)) stop("idVar should be a character! (eg. rankVar = 'rank')")
  
  colnames <- names(data)
  
  if (!(idVar %in% colnames)) stop("idVar not in dataframe!")
  if (!(rankVar %in% colnames)) stop("rankVar not in dataframe!")
  
  if (!is.null(covs.fix)) {
    if (!is.list(covs.fix)) stop("covs.fix should be a list!")
    if (any(!(unlist(covs.fix) %in% colnames))) stop("Some of covs.fix are missing in dataframe!")
  }
  
  if (!is.null(covs.het)) {
    if (!is.list(covs.het)) stop("covs.het should be a list!")
    if (any(!(unlist(covs.het) %in% colnames))) stop("Some of covs.het are missing in dataframe!")
  }
  
  if (!is.null(covsInt.fix)) {
    if (!is.list(covsInt.fix)) stop("covsInt.fix should be a list!")
    if (any(!(unlist(covsInt.fix) %in% colnames))) stop("Some of covsInt.fix are missing in dataframe!")
  }
  
  if (!is.null(covsInt.het)) {
    if (!is.list(covsInt.het)) stop("covsInt.het should be a list!")
    if (any(!(unlist(covsInt.het) %in% colnames))) stop("Some of covsInt.het are missing in dataframe!")
  }
  
  if (!is.null(FE)) {
    if (!is.character(FE)) stop("The option FE should be a character vector!")
    normalizeFE <- TRUE
    if (length(FE) > 1) normalizeFE <- TRUE
  }
}




