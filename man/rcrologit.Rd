% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rcrologit.R
\name{rcrologit}
\alias{rcrologit}
\title{Estimation and Inference for Random Coefficients Rank-Ordered Logit}
\usage{
rcrologit(
  dataPrep,
  Sigma = "diagonal",
  S = 50,
  approx.method = "MC",
  stdErr = "numerical",
  verbose = FALSE,
  control.opts = NULL
)
}
\arguments{
\item{dataPrep}{object of class 'rcrologit' prepared via \code{\link{dataPrep}}.}

\item{Sigma}{structure of the variance-covariance of the random coefficients. It must be one of
"diagonal" or "cholesky". Default is \code{Sigma="diagonal"}.}

\item{S}{integer denoting the number of simulations when approximating the integrals in the conditional
choice probabilities. Default is \code{S=50}.}

\item{approx.method}{string indicating the procedure to approximate the integrals in the conditional
choice probabilities}

\item{stdErr}{string denoting whether standard error should be estimated and how. Available options are
"numerical" (default) which uses numerical approximation of the Hessian matrix; "analytical" which uses
the closed form of the variance of the score and the hessian of the likelihood to estimate standard
errors robust to misspecification; "skip" which makes R skip the computation of standard errors}

\item{verbose}{if \code{TRUE} prints additional information in the console.}

\item{control.opts}{a list containing options to be passed to the underlying optimizer
\code{optim}.}
}
\value{
The function returns a list containing the following objects:
\item{b}{vector containing all the estimated parameters}
\item{bfix}{vector containing the estimated parameters corresponding to the \verb{fixed'' coefficients\} \\item\{bhet\}\{vector containing the estimated parameters corresponding to the }random'' coefficients}
\item{Lambda}{loading matrix of the shocks to the random coefficients}
\item{Sigma}{variance-covariance matrix of the estimated parameters}
}
\description{
Estimation and Inference for Random Coefficients Rank-Ordered Logit
}
\examples{
data <- rcrologit_data
dataprep <- dataPrep(data, idVar = "Worker_ID", rankVar = "rank",
                     altVar = "alternative",
                     covsInt.fix = list("Gender"),
                     covs.fix = list("log_Wage"), FE = c("Firm_ID"))
                     
rologitEst <- rcrologit(dataprep)

}
\seealso{
\code{\link{dataPrep}}
}
\author{
Chiara Motta, University of California Berkeley. \email{cmotta@berkeley.edu}

Filippo Palomba, Princeton University. \email{fpalomba@princeton.edu}
}