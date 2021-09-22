\name{Smooth}

\alias{Smooth}

\title{Smoothing data points}

\description{
Function for smoothing data points.
}

\usage{
Smooth(x,y,sigma,newx=NA,fac.extrap=1,method="quadratic")
}

\arguments{
  \item{x}{vector of x-coordinates of data points to be smoothed}
  \item{y}{vector of y-coordinates of data points to be smoothed (x and y must be of the same length)}
  \item{sigma}{standard deviation of Gaussian distribution used as weights for local quadratic or local linear regression}
  \item{newx}{optional vector of x-coordinates at which smoothed results and derivatives are to be calculated (if not specified, results are provided at the same locations as there is data available)}
  \item{fac.extrap}{calculate smoothed value only if data is available within \code{fac.extrap*sigma} (default is 1)}
  \item{method}{"quadratic" (default) indicates local quadratic regression, otherwise local linear regression is used}
}

\details{
Function for smoothing data points and estimating the derivative of the smoothed curve by local quadratic or optionally local linear regression. Local regression is implemented by using a Gaussian distribution of weights centered at the point at which the smoothed curve is evaluated. The smoothing parameter is the standard deviation of the Gaussian weights.
}

\value{
Returns a data frame containing:
  \item{x}{x-coordinates at which smoothed results and derivatives are available}
  \item{y}{smoothed results at the locations x}
  \item{ydot}{derivatives of smoothed results at the locations x}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}