\name{emulator.pri.var.normal}

\alias{emulator.pri.var.normal}

\title{Emulator: default model for prior variance}

\description{
Default model for prior variance (Gaussian stochastic process). Used as default for prior variance when constructing an emulator with \code{emulator.create}
}

\usage{
emulator.pri.var.normal(inp,sd,lambda,alpha=2,sd.smooth=0,rows=NA,cols=NA)}

\arguments{
  \item{inp}{input}
  \item{sd}{standard deviation}
  \item{lambda}{correlation length}
  \item{alpha}{exponent of distance to model decay of correlation}
  \item{sd.smooth}{smoothing}
  \item{rows}{rows (optional)}
  \item{cols}{columns (optional)}
}

\value{
Returns matrix of prior variance.
}

\seealso{
\code{\link{emulator.evaluate}}, \code{\link{emulator.create}}, \code{\link{emulator.pri.mean.lin}}, \code{\link{evaluate.emulator.GASP}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}