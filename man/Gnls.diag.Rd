\name{Gnls.diag}

\alias{Gnls.diag}

\title{Generalized nonlinear least squares regression diagnostics}

\description{
Calculate regression diagnostics for generalized nonlinear least squares.
}

\usage{
Gnls.diag(res.gnls,par.inc)
}

\arguments{
  \item{res.gnls}{output from \code{Gnls}}
  \item{par.inc}{increments of parameters used to approximate the derivatives}
}

\value{
Returns a list containing:
  \item{par.inc}{parameter increments used to calculate \code{V}}
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