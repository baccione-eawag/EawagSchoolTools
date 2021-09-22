\name{contourpdf}

\alias{contourpdf}

\title{Plotting functions (eawagSummerSchoolTools)}

\description{
Contour plot of a two-dimensional probability density function.
}

\usage{
contourpdf(calcpdf_mv,norm=T,xlim=c(-3,3),ylim=c(-3,3),
           levels=c(0.05,0.5,0.95),res=20,lty="solid",...)
}

\arguments{
  \item{calcpdf_mv}{function to calculate the log pdf values at a set of locations specified by its first argument}
  \item{norm}{TRUE if the probability density is normalized, FALSE if integration over the domain given by xlim and ylim should be used for normalization}
  \item{xlim}{bounds of the integration range for the first variable}
  \item{ylim}{bounds of the integration range for the second variable}
  \item{levels}{vector of probabilities to be contained in the contour line}
  \item{res}{resolution of grid used to countour (number of points in each dimension)}
  \item{lty}{line type of contour lines}
  \item{...}{additional arguments passed to \code{calcpdf_mv}}
}

\value{
Contour plot. Returns integral of the probability density over the given range at the given resolution.
}

\seealso{
\code{\link{Plot.res}}, \code{\link{Plot.vars}}, \code{\link{Plot.chains}}, \code{\link{Plot.margs}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}