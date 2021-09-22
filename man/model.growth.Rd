\name{model.growth}

\alias{model.growth}

\title{Example: Growth model}

\description{
Defines a simple growth model, which is used in the exercises for the Eawag Summer School in Environmental Systems Analysis.
}

\usage{
model.growth(par,L)
}

\arguments{
  \item{par}{named vector of parameters}
  \item{L}{observation layout}
}

\details{
Calculates the growth of micro-organisms on a substrate in a batch reactor:

dC_M/dt = mu * C_S/(K+C_S) * C_M - b*C_M

dC_S/dt = -(mu/Y) * C_S/(K+C_S) * C_M

State variables are:

\tabular{ll}{
  \code{C_M} \tab concentration of micro-organisms in the reactor \cr
  \code{C_S} \tab concentration of substrate in the reactor \cr
}

Parameters are:

\tabular{ll}{
  \code{mu} \tab maximum growth rate of micro-organisms \cr
  \code{K} \tab half-concentration of growth rate with respect to substrate \cr
  \code{b} \tab rate of death and respiration processes of micro-organisms \cr
  \code{Y} \tab yield of growth process \cr
  \code{C_M_ini} \tab initial concentration of micro-organisms \cr
  \code{C_S_ini} \tab initial concentration of substrate \cr
}

The observation layout \code{L} contains the variable names and time points at which the model is to be evaluated. \code{L} should be a data frame with two columns. The first column contains the names of the variables to be evaluated, and the second column contains the time points (or other evaluation steps) for the model. See function \code{checkLayout} in package \code{bayesBias} for more details. For example, for the growth model:

\tabular{cc}{
  var \tab t \cr
  C_M \tab 1 \cr
  C_M \tab 2 \cr
  C_M \tab 3 \cr
  ... \tab ... \cr
  C_S \tab 1 \cr
  C_S \tab 2 \cr
  C_S \tab 3 \cr
  ... \tab ... \cr
}
}

\value{
Returns a named vector of model output (concentrations). Names consist of the variable name and evaluation point (time), separated by an underscore (for example \code{C_M_1}).
}

\seealso{
\code{\link{model.monod}}, \code{\link{model.monod.external}}, and function \code{checkLayout} in package \code{bayesBias}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>, Lukas M. Weber
}