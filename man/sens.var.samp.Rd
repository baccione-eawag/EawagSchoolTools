\name{sens.var.samp}

\alias{sens.var.samp}

\title{Variance-based sensitivities of a deterministic model}

\description{
Function to calculate the first order variance-based sensitivity coefficients from a parameter sample and the corresponding result sample.
}

\usage{
sens.var.samp(parsamp,ressamp,nbin=NA,method="smooth",order=2,
              bandwidth=1,span=0.75,sd.rel=0.1,plot=F)
}

\arguments{
  \item{parsamp}{matrix containing the parameter sample (each row corresponds to a sample point)}
  \item{ressamp}{vector (of only one result per parameter sample point) or matrix of results corresponding to the parameter sample (each row provides the results corresponding to the parameter values in the same row of parsamp)}
  \item{nbin}{number of quantile intervals for which the conditional means are calculated, default is the square root of the sample size}
  \item{method}{"smooth", "loess", "glkerns", "lokerns", "lpepa", "lpridge": routine to be used for smoothing}
  \item{order}{order of local regression polynomial or of kernel}
  \item{bandwidth}{\code{method="lpepa"} or \code{method="lpdidge"} only: bandwidth of the smoothing algorithm}
  \item{span}{\code{method="loess"} only: fraction of points used for local regression}
  \item{sd.rel}{\code{method="smooth"} only: standard deviation of normal distribution of smoothing algorithm relative to the 99\% quantile interval}
  \item{plot}{logical variable indicating if a scatter plot of the relationship between parameter and model output should be plotted}
}

\value{
Returns a list with two elements:
  \item{var}{vector with the total variance of each column (each model output) of the ressamp matrix}
  \item{var.cond.E}{matrix with the variance of the conditional expected value of each parameter (columns) and each model output (rows)}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}