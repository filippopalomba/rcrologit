################################################################################
#' Print Method for Random Coefficients Rank-Ordered Logit
#'
#' @description The print method for for random coefficients rank-ordered logit fitted objects.
#'
#' @param x Class "rcrologit" object, obtained by calling \code{\link{rcrologit}}.
#' @param printFE Whether fixed effects coefficients and standard errors should be displayed.
#' @param ... Other arguments.
#'
#' @return No return value, called to print \code{\link{rcrologit}} results.
#'
#' @author
#' Chiara Motta, University of California Berkeley. \email{cmotta@berkeley.edu}
#'
#' Filippo Palomba, Princeton University. \email{fpalomba@princeton.edu}
#'
#' @seealso \code{\link{rcrologit}} for estimation and inference for random coefficients rank-ordered logit.
#'
#' Supported methods: \code{\link{print.rcrologit}}.
#'
#' @export
#'

print.rcrologit <- function(x, printFE = FALSE, ...) {

  rCoefs <- x$param.spec$model == 'rcrologit'
  K_il.fix <- x$param.spec$K_il.fix
  K_i.fix <- x$param.spec$K_i.fix
  K.fix <- x$param.spec$K.fix
  Kfe <- x$param.spec$Kfe
  K_il.het <- x$param.spec$K_il.het
  K_i.het <- x$param.spec$K_i.het
  
  cat("\n")
  if (rCoefs == FALSE) {
    cat(paste0("Rank Ordered Logit - Results\n"))
  } else {
    cat(paste0("Random Coefficients Rank Ordered Logit - Results\n"))
    K.lambda <- nrow(x$Lambda) * (nrow(x$Lambda) + 1) / 2
  }
  cat("\n")

  # print (if present) fixed coefficients
  if (K_il.fix > 0 | K_i.fix > 0) {
    K.fixx <- K.fix
    if (printFE == FALSE) K.fixx <- K.fix - Kfe 
      
    bFixprint <- as.matrix(x$bfix[1:(K.fixx)])
    seFixprint <- as.matrix(sqrt(diag(x$Sigma)[1:(K.fixx)]))
    namesFix <- names(x$bfix[1:(K.fixx)])
    lb <- bFixprint - qnorm(0.975) * seFixprint
    ub <- bFixprint + qnorm(0.975) * seFixprint
    
    FixPrint <- round(cbind(bFixprint, seFixprint, lb, ub), 3)
    rownames(FixPrint) <- namesFix
    colnames(FixPrint) <- c("Coef.", "Std. Err.", "Lower Bound", "Upper Bound")
    
    cat("Fixed Coefficients:\n")
    print(FixPrint, col.names = FALSE)
    
  }

  # print (if present) random coefficients
  if (K_il.het > 0 | K_i.het > 0) {
    
    # prepare mean coefficients 
    fixParams <- K.fix + 1
    
    bRanPrint <- as.matrix(x$bhet)
    seRanPrint <- as.matrix(sqrt(diag(x$Sigma)[fixParams:(fixParams + nrow(bRanPrint) - 1)]))
    namesHet <- names(x$bhet)

    lb <- bRanPrint - qnorm(0.975) * seRanPrint
    ub <- bRanPrint + qnorm(0.975) * seRanPrint
    
    RanPrint <- round(cbind(bRanPrint, seRanPrint, lb, ub), 3)
    rownames(RanPrint) <- namesHet
    colnames(RanPrint) <- c("Coef.", "Std. Err.", "Lower Bound", "Upper Bound")
    
    # prepare standard error coefficients
    SigmaB <- x$Lambda %*% t(x$Lambda)

    if (x$param.spec$Sigma == "diagonal") {
      SigmaBvech <- diag(SigmaB)
      SigmaN <- strsplit(rownames(x$Lambda),"\\.")[[1]]
      SigmaN <- unlist(purrr::map(strsplit(rownames(x$Lambda), "\\."), 2))
    } else {
      k <- nrow(SigmaB)
      SigmaRN <- rownames(x$Lambda)
      SigmaCN <- colnames(x$Lambda)
      SigmaBvech <- c()
      SigmaN <- c()
      for (i in seq_len(k)) {
        for (j in seq_len(i)) {
          SigmaBvech <- c(SigmaBvech, SigmaB[i, j])
          if (i == j) {
            name <- strsplit(SigmaRN[i],"\\.")[[1]][2]
          } else {
            name <- paste0(strsplit(SigmaRN[i],"\\.")[[1]][2],"-",
                           strsplit(SigmaRN[j],"\\.")[[1]][2])
          }
          SigmaN <- c(SigmaN, name)
        }
      }
    }

    names(SigmaBvech) <- SigmaN
    seRanSPrint <- as.matrix(sqrt(diag(x$Sigma)[(fixParams + nrow(bRanPrint)):ncol(x$Sigma)]))

    lb <- SigmaBvech - qnorm(0.975) * seRanSPrint
    ub <- SigmaBvech + qnorm(0.975) * seRanSPrint

    RanSPrint <- round(cbind(SigmaBvech, seRanSPrint, lb, ub), 3)
    rownames(RanSPrint) <- SigmaN
    colnames(RanSPrint) <- c("Coef.", "Std. Err.", "Lower Bound", "Upper Bound")
    
    cat("\n")
    cat("Random Coefficients - Mean:\n")
    print(RanPrint, col.names = FALSE)

    cat("\n")
    cat("Random Coefficients - Variances:\n")
    print(RanSPrint, col.names = FALSE)
    
  }
  
}

################################################################################
#' Summary Method for Random Coefficients Rank-Ordered Logit
#'
#' @description The summary method for for random coefficients rank-ordered logit fitted objects.
#'
#' @param object Class "rcrologit" object, obtained by calling \code{\link{rcrologit}}.
#' @param printFE Whether fixed effects coefficients and standard errors should be displayed.
#' @param ... Other arguments.
#'
#' @return No return value, called to print \code{\link{rcrologit}} results.
#'
#' @author
#' Chiara Motta, University of California Berkeley. \email{cmotta@berkeley.edu}
#'
#' Filippo Palomba, Princeton University. \email{fpalomba@princeton.edu}
#'
#' @seealso \code{\link{rcrologit}} for estimation and inference for random coefficients rank-ordered logit.
#'
#' Supported methods: \code{\link{print.rcrologit}}.
#'
#' @export
#'

summary.rcrologit <- function(object, printFE = FALSE, ...) {
  
  cat("\n")
  if (object$param.spec$model == 'rcrologit') {
    cat(paste0("Rank Ordered Logit - Setup\n"))
  } else {
    cat(paste0("Random Coefficients Rank Ordered Logit - Setup\n"))
  }
  cat("\n")
  
  K_il.fix <- object$param.spec$K_il.fix                     # number of unit-alternative varying covs with fixed taste
  K_il.het <- object$param.spec$K_il.het                     # number of unit-alternative varying covs with het taste
  K_i.fix <- object$param.spec$K_i.fix                       # number of unit varying covs with fixed taste
  K_i.het <- object$param.spec$K_i.het                       # number of unit varying covs with het taste
  J <- object$param.spec$J                                   # number of alternatives 
  N <- object$param.spec$N                                   # number of observations  
  if (object$param.spec$Kfe > 0) {
    Kfestr <- "Yes"
  } else {
    Kfestr <- "No"
  }
  cat("--------------------------------------------------------------------\n")
  cat(paste("Number of units:                                           ", N,"\n", sep = ""))
  cat(paste("Number of alternatives:                                    ", J,"\n", sep = ""))
  cat(paste("Unit-alternative varying covariates with fixed taste:      ", K_il.fix, "\n", sep = ""))
  cat(paste("Unit varying covariates with fixed taste:                  ", K_i.fix,"\n", sep = ""))
  cat(paste("Unit-alternative varying covariates with random taste:     ", K_il.het, "\n", sep = ""))
  cat(paste("Unit varying covariates with random taste:                 ", K_i.het,"\n", sep = ""))
  cat(paste("Unit fixed effects:                                        ", Kfestr,"\n", sep = ""))
  
  cat("--------------------------------------------------------------------\n")
  
  print.rcrologit(object)
}
