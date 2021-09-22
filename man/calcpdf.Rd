\name{calcpdf}

\alias{calcpdf}

\title{Densities of univariate distributions}

\description{
Convenience function to return probability densities of several univariate distributions.
}

\usage{
calcpdf(x,distpar,log=FALSE)
}

\arguments{
  \item{x}{vector at which probability densities will be evaluated}
  \item{distpar}{distribution type and parameters (see details)}
  \item{log}{if TRUE returns log of density, otherwise density}
}

\details{
Argument \code{distpar} should be a character vector containing the distribution name (first element in the vector) and parameter values (subsequent elements; these are converted to numeric values automatically).

Implemented distribution types and their required parameters are:

"Uniform" : min, max \cr
"Normal" : mean, sd \cr
"NormalTrunc" (truncated normal distribution) : mean, sd, min, max \cr
"Lognormal" : mean, sd (note R parameters are mean and sd of the log of the random variable) \cr
"LognormalTrunc" : mean, sd, min, max (note R parameters are mean and sd of the log of the random variable) \cr
"Inv" (inverse distribution) : min, max \cr
"Exponential" : mean \cr
}

\value{
Probability density (or log density) at the sample points in x.
}

\seealso{
\code{\link{calcpdf_mv}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}