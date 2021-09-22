\name{resid.diag.boxcox}

\alias{resid.diag.boxcox}

\title{Plotting functions (eawagSummerSchoolTools)}

\description{
Plot residual diagnostic plots for given Box-Cox transformation parameters.
}

\usage{
resid.diag.boxcox(obs,calc,lambda1=1,lambda2=1)
}

\arguments{
  \item{obs}{vector of observations}
  \item{calc}{vector of calculated results}
  \item{lambda1}{value of parameter \code{lambda1} (exponent) of Box-Cox transformation}
  \item{lambda2}{value of parameter \code{lambda2} (offset) of Box-Cox transformation}
}

\value{
Diagnostic plots.
}

\seealso{
\code{\link{resid.diag}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}