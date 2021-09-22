\name{intpol.irregular}

\alias{intpol.irregular}

\title{Interpolation on an irregular grid}

\description{
Function to perform interpolation on an irregular grid.
}

\usage{
intpol.irregular(x,x.data,y.data,method=c("neighbor","invdist"),
                 scales.x=NA,fact=1)
}

\arguments{
  \item{x}{vector of values of independent variable at which to interpolate; \code{xout} can also be a matrix, interpolation is then for each row}
  \item{x.data}{matrix of independent independent variables; each row of \code{x.data} belongs to one input specification (\code{x.data} can be a vector if there is only one input dimension)}
  \item{y.data}{vector of values of the dependent variable; one value for each row of \code{x.data}}
  \item{method}{"neighbor" does nearest neighbor interpolation (piecewise constant function of value at nearest data point); "invdist" weights the provided dependent values according to the inverse of the distance to \code{xout} squared}
  \item{scales.x}{values of \code{x} and of \code{x.data} are scaled by these scales before interpolation; default is all scales equal to 1}
  \item{fact}{scaling factor}
}

\value{
Returns interpolated values.
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}