\name{Impsamp}

\alias{Impsamp}

\title{Importance sampling}

\description{
Calculate an importance sample of a distribution.
}

\usage{
Impsamp(log.pdf,z,z.log.pdf,...)
}

\arguments{
  \item{log.pdf}{log probability density function}
  \item{z}{parameter samples}
  \item{z.log.pdf}{log probability density at the sample points}
  \item{...}{additional arguments passed to \code{log.pdf}}
}

\value{
Returns a list containing:
  \item{z}{values of z}
  \item{w}{values of w}
  \item{ess}{values of ess}
  \item{log.pdf.dist}{values of log probability density}
  \item{log.pdf.samp}{log probability density at the sample points (\code{z.log.pdf})}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}