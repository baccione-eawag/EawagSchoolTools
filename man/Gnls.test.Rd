\name{Gnls.test}

\alias{Gnls.test}

\title{Generalized nonlinear least squares test statistics}

\description{
Calculate test statistics for generalized nonlinear least squares.
}

\usage{
Gnls.test(res.gnls,V,par,y.det)
}

\arguments{
  \item{res.gnls}{output from Gnls}
  \item{V}{matrix of partial derivatives of model results with respect to parameters at the estimated parameter values}
  \item{par}{named parameter vector with parameter estimates}
  \item{y.det}{deterministic model results corresponding to the estimated parameters}
}

\value{
Returns a list containing:
  \item{V}{matrix of partial derivatives of model results with respect to parameters at the estimated parameter values}
  \item{sd.par}{standard errors of the parameters assuming known error variances}
  \item{var.par}{estimated variance-covariance matrix of the parameters assuming known error variances}
  \item{corr.par}{estimated correlation matrix of the parameters}
  \item{sd.rel}{estimated correction factor in standard deviations of the error model}
  \item{sd.par.rel}{standard errors of the parameters when estimating a common factor in error variances}
  \item{var.par.rel}{estimated variance-covariance matrix of the parameters when estimating a common factor in error variances}
}

\seealso{
\code{\link{Gnls}}, \code{\link{Gnls.diag}}, \code{\link{Gnls.test}}, \code{\link{Gnls.predict}}, \code{\link{Confint}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}