\name{emulator.pri.mean.lin}

\alias{emulator.pri.mean.lin}

\title{Emulator: default model for prior mean}

\description{
Default model for prior mean (linear model). Used as default for prior mean when constructing an emulator with \code{emulator.create}.
}

\usage{
emulator.pri.mean.lin(inp,par)
}

\arguments{
  \item{inp}{input matrix}
  \item{par}{parameter vector}
}

\value{
Returns vector of prior mean.
}

\seealso{
\code{\link{emulator.evaluate}}, \code{\link{emulator.create}}, \code{\link{emulator.pri.var.normal}}, \code{\link{evaluate.emulator.GASP}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}