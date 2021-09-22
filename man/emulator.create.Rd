\name{emulator.create}

\alias{emulator.create}

\title{Create emulator}

\description{
Function to construct an emulator based on design data and parameters.
}

\usage{
emulator.create(inp.design,res.design,par,sd,lambda,alpha=2,sd.smooth=0,
                pri.mean=emulator.pri.mean.lin,pri.var=emulator.pri.var.normal)
}

\arguments{
  \item{inp.design}{input design}
  \item{res.design}{result design}
  \item{par}{parameter}
  \item{sd}{parameter}
  \item{lambda}{parameter}
  \item{alpha}{parameter}
  \item{sd.smooth}{parameter}
  \item{pri.mean}{function for prior mean; default \code{emulator.pri.mean.lin}}
  \item{pri.var}{function for prior variance; default \code{emulator.pri.var.normal}}
}

\value{
Returns an emulator as a list with class label "emulator.GASP".
}

\seealso{
\code{\link{emulator.evaluate}}, \code{\link{emulator.pri.mean.lin}}, \code{\link{emulator.pri.var.normal}}, \code{\link{evaluate.emulator.GASP}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}