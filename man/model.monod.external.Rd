\name{model.monod.external}

\alias{model.monod.external}

\title{Example: Monod model (external)}

\description{
Defines a simple Monod model, which is used in the exercises for the Eawag Summer School in Environmental Systems Analysis. The "external" version \code{model.monod.external} calls an external program (\code{model_monod.exe}) to calculate the model results, while the "internal" version \code{model.monod} implements the model directly in R.
}

\usage{
model.monod.external(par,L)
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

The "external" version of the model \code{model.monod.external} implements the same equation as the "internal" version \code{model.monod}. The external version writes parameter values and inputs to files, then calls an external program (in this case \code{model_monod.exe}) that reads the parameter and input files, calculates the results, and writes the results to a text file. The external model then reads the results from the text file, and returns the values.

The external model \code{model.monod.external} can be used as a template to implement calls to other external simulation programs.
}

\value{
Returns a named vector of model output (growth rates). Names consist of the variable name and evaluation point (concentration), separated by an underscore (for example \code{r_0.1}).
}

\seealso{
\code{\link{model.monod}}, \code{\link{model.growth}}, and function \code{checkLayout} in package \code{bayesBias}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>, Lukas M. Weber
}