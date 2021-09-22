\name{calcpdf_mv}

\alias{calcpdf_mv}

\title{Log density of multivariate normal, lognormal, or product of 1d distributions}

\description{
Calculate the logarithm of the probability density function of a multivariate normal or lognormal distribution, or of a product of independent one-dimensional distributions.
}

\usage{
calcpdf_mv(z,dist="normal",mean=0,sd=1,cor=0,cor.inv=NA,log=TRUE,
           distdef=NA,file=NA)
}

\arguments{
  \item{z}{vector, matrix or data frame at which the logarithm of the probability density function will be evaluated}
  \item{dist}{distribution type: "Normal", "Lognormal", or "Indep"}
  \item{mean}{vector of means}
  \item{sd}{vector of standard deviations}
  \item{cor}{correlation matrix of the distribution}
  \item{cor.inv}{inverse of correlation matrix of the distribution (alternative input to cor, saves computation time for repeated calls)}
  \item{log}{if TRUE returns log of pdf, otherwise pdf}
  \item{distdef}{distribution definition for independent one-dimensional distributions}
  \item{file}{optional file to save results}
}

\value{
Probability density (or log density) at the sample points in z.
}

\seealso{
\code{\link{calcpdf}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}