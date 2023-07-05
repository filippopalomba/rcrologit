##########################################################################
##########################################################################
# computes the likelihood of the rank-ordered logit
loglkld <- function(b0, X, J) {

  eXb <- lapply(X, function(x) exp(x%*%b0))
  
  lkld <- 0
  for (j in seq_len(J-1)) {
    lkld <- lkld + log(eXb[[j]] / unlist(Reduce(`+`, eXb[c(j:J)])))
  }  
    
  return(-sum(lkld))
}


##########################################################################
# computes the likelihood of the random coefficients rank-ordered logit
loglkldRC <- function(b0, X, J, K.fix, K.het.mu, K.het.lam, Sigma, bias.corr,
                      approx.method, S, epsMC, Ncores) {

  # unpack b0 to pass it to mvnorm
  k1 <- K.fix                            
  k2 <- K.fix + 1
  k3 <- k1 + K.het.mu
  k4 <- k3 + 1
  k5 <- k3 + K.het.lam
  bfix <- b0[1:k1]                                # first sub-component is b for fixed taste covs (including FE)
  bhet <- b0[k2:k3]                               # second: mean het taste 
  bLam <- vech2mat(b0[k4:k5], K.het.mu, Sigma)    # third: varcov of shocks het taste 
  
  ccp <- ccpGet(bfix, bhet, bLam, X, J, bias.corr, approx.method, S, epsMC, Ncores) 
  lkld <- -sum(ccp)
  
  return(lkld)
}

##########################################################################
# estimates the conditional choice probabilities
ccpGet <- function(bfix, bhet, bLam, X, J, bias.corr, approx.method, S, epsMC, Ncores) {
  
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
    if (Ncores == 1) {
      ccp <- apply(b_i, 1, function(x) ccpROLogit(X, x, J))
    } else {
      b_list <- lapply(seq_len(nrow(b_i)), function(i) b_i[i,])
      tmp <- parallel::mclapply(b_list, function(x) ccpROLogit(X, x, J), 
                                mc.cores=Ncores)
      ccp <- matrix(unlist(tmp), nrow=nrow(X[[1]]), ncol=S, byrow=FALSE)
    }
    
    if (is.null(dim(ccp))) ccp <- t(as.matrix(ccp)) # handles SE computation
    
    fhat <- rowMeans(ccp)
    
    if (bias.corr == TRUE) {

      fhatVar <- rowMeans((ccp - fhat)^2)
      bc <- 0.5 * fhatVar / fhat^2 
      
    } else {
      
      bc <- 0
      
    }
    
    ccp <- log(fhat) + bc
  }
  
  return(ccp)
}

##########################################################################
# compute rank ordered logit ccp
ccpROLogit <- function(X, b, J) {
  eXb <- lapply(X, function(x) exp(x%*%as.matrix(b)))
  rol <- 1
  for (j in seq_len(J-1)) {
    rol <- rol * (eXb[[j]] / unlist(Reduce(`+`, eXb[c(j:J)])))
  }
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
seGet <- function(Xlist, b, rCoefs, robust, pars = NULL) {

  ###########################
  # Standard rologit  
  ###########################
  if (rCoefs == FALSE) {
    
    # analytical formula to get robust (or non-robust) standard errors
    
    N <- pars$N
    J <- pars$J

    # Jacobian
    eXb <- lapply(Xlist, function(x) exp(x%*%as.matrix(b)))
    XeXb <- lapply(c(1:J), function(i) Xlist[[i]] * c(eXb[[i]]))
    
    Jac <- 0
    for (j in seq_len(J-1)) {
      Jac <- Jac + Xlist[[j]] - unlist(Reduce(`+`, XeXb[c(j:J)])) / c(unlist(Reduce(`+`, eXb[c(j:J)])))
    }
    Jac <- var(Jac)
    
    # Hessian
    
    H <- 0
    for (i in seq_len(N)) {
      xlist <- lapply(Xlist, function(x) t(x[i,,drop=FALSE]))            # X vector for each alternative
      exblist <- lapply(xlist, function(x) exp(sum(x*b)))                # exp(X'b) for each alternative
      Xexblist <- lapply(c(1:J), function(j) xlist[[j]] * exblist[[j]])  # X * exp(X'b) for each alternative
      XXexblist <- lapply(c(1:J), function(j) xlist[[j]] %*% t(xlist[[j]]) * exblist[[j]])  # X * X' * exp(X'b) for each alternative
      
      Hes_i <- 0 
      for (j in seq_len(J-1)) {
        h1 <- unlist(Reduce(`+`, XXexblist[c(j:J)]))
        h2 <- sum(unlist(exblist[j:J]))
        den <- h2^2
        Xexb <- unlist(Reduce(`+`, Xexblist[c(j:J)]))
        h3 <- Xexb %*% t(Xexb)
        Hes_i <- Hes_i + (1/den) * (h1*h2 - h3)
      }
      
      H <- H + Hes_i / N
    }
      
    if (robust == TRUE) {
      
      Sigma <- solve(H) %*% Jac %*% solve(H) / N
      
    } else if (robust == FALSE) {
      
      Sigma <- solve(Jac) / N
      
    }
  }

  #############################
  # Random Coefficients rologit  
  #############################

  if (rCoefs == TRUE) {

    H <- numDeriv::hessian(loglkldRC, x = b, X = Xlist, J=pars$J, K.fix=pars$K.fix, K.het.mu=pars$K.het.mu,
                           K.het.lam=pars$K.het.lam, Sigma=pars$Sigma, bias.corr=pars$bias.corr, 
                           approx.method=pars$approx.method, S=pars$S, epsMC=pars$epsMC, Ncores=pars$Ncores)    
    
    # numerical approximation of hessian of likelihood function
    if (robust == FALSE) {

      Sigma <- solve(H)
      Jac <- NULL
      
    } else {
      # numerical approximation of the variance of the jacobian of likelihood function
      
      if(pars$verbose && pars$Ncores == 1) {
        cat("This might take several minutes! If you are not already doing it, consider using multiple cores!")
      }
      
      ccp <- ccpGet(bfix=pars$bfix, bhet=pars$bhet, bLam=pars$bLam, X=Xlist, J=pars$J, 
                    bias.corr=pars$bias.corr, approx.method=pars$approx.method, S=pars$S,
                    epsMC=pars$epsMC, Ncores=pars$Ncores)
            
      Jlist <- parallel::mclapply(c(1:pars$N), function(i) jacobianGet(i=i, Xlist=Xlist, b=b, pars=pars),
                                  mc.cores = pars$Ncores)
      Jmat <- Reduce(rbind, Jlist)
      Jac <- var(Jmat * sqrt(pars$N)) # rescale by sqrt(N) 
      
      Sigma <- solve(H) %*% Jac %*% solve(H) 
    }
  }

  return(list(Sigma = Sigma, Jac = Jac, H = H)) 
}


jacobianGet <- function(i, Xlist, b, pars) {
  
  xdata <- lapply(Xlist, function(x) x[i, , drop=F])
  J <- numDeriv::jacobian(loglkldRC, x = b, X = xdata, J=pars$J, K.fix=pars$K.fix, K.het.mu=pars$K.het.mu,
                          K.het.lam=pars$K.het.lam, Sigma=pars$Sigma, bias.corr=pars$bias.corr, 
                          approx.method=pars$approx.method, S=pars$S, epsMC=pars$epsMC, Ncores=1)

  return(J)
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

}




