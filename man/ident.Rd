\name{ident}

\alias{ident}

\title{Identifiability measures}

\description{
Function to calculate identifiability measures.
}

\usage{
ident(sen,delta.par=0,scale=0,max.subset.size=0)
}

\arguments{
  \item{sen}{matrix of model sensitivities (partial derivatives of model outcomes with respect to model parameters: dy/dpar); the columns of sen refer to different model parameters, the rows to different model outcomes}
  \item{delta.par}{model parameter uncertainty ranges (length equal to the number of columns of sen); if zero, then all ranges are assumed to be unity}
  \item{scale}{scaling factors of model results (if zero, then all scaling factors are assumed to be unity)}
  \item{max.subset.size}{parameter for maximum subset size}
}

\details{
This function calculates a parameter sensitivity ranking and collinearity indices for a series of parameter combinations based on linear sensitivity functions of a model, parameter uncertainty ranges and scale factors of model results.

The function is a simplified version of the program "ident" available at http://www.ident.eawag.ch.
}

\value{
Returns a list containing: \code{delta.msqr}, \code{collind}.
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}