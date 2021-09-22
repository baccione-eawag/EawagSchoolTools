\name{evaluate.emulator.GASP}

\alias{evaluate.emulator.GASP}

\title{Class method to evaluate emulator}

\description{
Class method to evaluate an emulator with class label "emulator.GASP".
}

\usage{
evaluate.emulator.GASP(emulator,inp)
}

\arguments{
  \item{emulator}{an emulator that has been constructed using \code{emulator.create}}
  \item{inp}{input}
}

\value{
Returns a list containing:
  \item{inp}{input}
  \item{y}{output}
  \item{var}{variance}
}

\seealso{
\code{\link{emulator.evaluate}}, \code{\link{emulator.create}}, \code{\link{emulator.pri.mean.lin}}, \code{\link{emulator.pri.var.normal}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}