\name{Smooth.fun}

\alias{Smooth.fun}

\title{Smoothing piecewise linear functions}

\description{
Function for smoothing piecewise linear functions.
}

\usage{
Smooth.fun(data,z,sigma,newz=NA,newx=NA,fac.extrap=1,method="quadratic")
}

\arguments{
  \item{data}{list of matrices specifying piecewise linear functions: for each function, the independent variable x must be provided in the first column, the dependent variable y in the second}
  \item{z}{vector of z-coordinates corresponding to the functions}
  \item{sigma}{standard deviation of Gaussian distribution used as weights for local quadratic regression in z-coordinates (see function \code{Smooth})}
  \item{newz}{optional vector of z-coordinates at which smoothed functions are to be calculated (if not specified, results are provided at the same locations as there is data available)}
  \item{newx}{optional vector of x-coordinates at which smoothed function values are to be calculated}
  \item{fac.extrap}{calculate smoothed value only if data is available within \code{fac.extrap*sigma} (default is 1)}
  \item{method}{"quadratic" (default) indicates local quadratic regression, otherwise local linear regression is used}
}

\details{
The smoothing parameter is the standard deviation of the Gaussian weights.
}

\value{
Returns a list of data frames with columns:
  \item{x}{x-coordinates at which smoothed results are available}
  \item{y}{smoothed results at the locations x}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}