\name{randsamp}

\alias{randsamp}

\title{Random sampling from multivariate normal, lognormal, or product of 1d distributions}

\description{
Generate a random sample from a multivariate normal or lognormal distribution, or from a product of independent one-dimensional distributions.
}

\usage{
randsamp(sampsize=1,dist="Normal",mean=0,sd=1,cor=NA,distdef=NA,file=NA)
}

\arguments{
  \item{sampsize}{sample size}
  \item{dist}{distribution type: "Normal", "Lognormal", or "Indep"}
  \item{mean}{vector of means}
  \item{sd}{vector of standard deviations}
  \item{cor}{correlation matrix of the distribution}
  \item{distdef}{distribution definition for independent one-dimensional distributions}
  \item{file}{optional file to save results}
}

\value{
Returns a list containing:
  \item{mean}{vector of means}
  \item{sd}{vector of standard deviations}
  \item{corr}{correlation matrix of the distribution}
  \item{sampsize}{sample size}
  \item{sample}{data frame of parameter samples (each row corresponds to a draw)}
  \item{logpdf}{vector of values of log probability density at the sample points}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}