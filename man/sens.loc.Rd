\name{sens.loc}

\alias{sens.loc}
\alias{sens.loc.explicitpars}

\title{Local sensitivities of a deterministic model}

\description{
Calculate local sensitivity matrix of a deterministic model (matrix of partial derivatives of model results with respect to parameters), based on a simple discretization of the derivatives.

\code{sens.loc.explicitpars} is a modified version of \code{sens.loc} which calls a model with explicitly listed arguments (this is primarily for internal use when calling \code{nls}).
}

\usage{
sens.loc(par,model,par.inc=0.01*par,...)

sens.loc.explicitpars(par,model.name,par.inc=0.01*par,x)
}

\arguments{
  \item{par}{named parameter vector}
  \item{model}{function representing the model; typically returns a vector, but could also return a matrix}
  \item{model.name}{function with explicit parameter arguments when using \code{sens.loc.explicitpars}}
  \item{par.inc}{increments of parameters used to approximate the derivatives}
  \item{x}{additional parameter arguments when using \code{sens.loc.explicitpars}}
  \item{...}{additional arguments passed to \code{model}}
}

\value{
Matrix of partial derivatives of model results (rows) with respect to parameters (columns), or 3-dimensional array of partial derivatives if model output is a matrix.
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}