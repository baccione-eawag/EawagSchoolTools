\name{emulator.evaluate}

\alias{emulator.evaluate}

\title{Evaluate emulator}

\description{
Function to evaluate an emulator at new points. An emulator can be constructed using \code{emulator.create}.
}

\usage{
emulator.evaluate(x,...)
}

\arguments{
  \item{x}{an emulator that has been constructed using \code{emulator.create}}
  \item{...}{new points passed to the emulator \code{x}}
}

\value{
Output of the emulator at the new points.
}

\seealso{
\code{\link{emulator.create}}, \code{\link{emulator.pri.mean.lin}}, \code{\link{emulator.pri.var.normal}}, \code{\link{evaluate.emulator.GASP}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}