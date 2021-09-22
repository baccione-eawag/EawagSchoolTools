\name{model.monod}

\alias{model.monod}

\title{Example: Monod model}

\description{
Defines a simple Monod model, which is used in the exercises for the Eawag Summer School in Environmental Systems Analysis.
}

\usage{
model.monod(par,L)
}

\arguments{
  \item{par}{named vector of parameters}
  \item{L}{observation layout}
}

\details{
The Monod model calculates the growth rate as a function of substrate concentration:

r = r_max*C / (K + C)

Parameters are:

\tabular{ll}{
  \code{r_max} \tab maximum growth rate \cr
  \code{K} \tab half-saturation concentration \cr
}

The observation layout \code{L} contains the concentrations C at which the growth rate is to be evaluated. \code{L} should be a data frame with two columns. The first column contains the names of the variables to be evaluated, and the second column contains the time points (or other evaluation steps) for the model. See function \code{checkLayout} in package \code{bayesBias} for more details. For example, for the Monod model:

\tabular{cc}{
  var \tab C \cr
  r   \tab 0.1 \cr
  r   \tab 0.2 \cr
  r   \tab 0.3 \cr
  ... \tab ... \cr
}
}

\value{
Returns a named vector of model output (growth rates). Names consist of the variable name and evaluation point (concentration), separated by an underscore (for example \code{r_0.1}).
}

\seealso{
\code{\link{model.monod.external}}, \code{\link{model.growth}}, and function \code{checkLayout} in package \code{bayesBias}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>, Lukas M. Weber
}