\name{Gnls.predict}

\alias{Gnls.predict}

\title{Generalized nonlinear least squares predictions}

\description{
Calculate predictions based on generalized nonlinear least squares regression results.
}

\usage{
Gnls.predict(diag.gnls,newx,level=0.95,var=NA)
}

\arguments{
  \item{diag.gnls}{output from \code{Gnls.diag}}
  \item{newx}{new model input x at which confidence intervals are to be calculated}
  \item{level}{probability level of confidence intervals (default 0.95)}
  \item{var}{error variance-covariance matrix at new model input (only variances used)}
}

\value{
Returns a list containing:
  \item{x}{x values at which model results are calculated}
  \item{y.pred}{predicted model results}
  \item{confint}{confidence intervals of deterministic model results}
  \item{confint.rel}{confidence intervals of deterministic model results when estimating a common factor in error variances}
  \item{predint}{prediction intervals (only if error variance-covariance matrix is provided)}
  \item{predint.rel}{prediction intervals (only if error variance-covariance matrix is provided) when estimating a common factor in error variances}
}

\seealso{
\code{\link{Gnls}}, \code{\link{Gnls.diag}}, \code{\link{Gnls.test}}, \code{\link{Gnls.predict}}, \code{\link{Confint}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}