\name{Loglikeli}

\alias{Loglikeli}

\title{Log likelihood}

\description{
Function implementing a simple standard log likelihood function for a model provided by a deterministic function. This function provides the density of a multivariate normal distribution centered at the results of the deterministic model.

This likelihood implementation serves as a template for implementing more specific likelihood functions. It is called by the function \code{Logposterior}.
}

\usage{
Loglikeli(par,model,L,y,error.sd=NA,error.cor.inv=NA,...)
}

\arguments{
  \item{par}{parameters}
  \item{model}{model function}
  \item{L}{observation layout}
  \item{y}{data}
  \item{error.sd}{used for variance-covariance structure}
  \item{error.cor.inv}{used for variance-covariance structure}
  \item{...}{additional arguments passed to \code{model}}
}

\details{
If there exists a component of the parameter vector labeled \code{sd_rel} instead of the provided standard deviations, \code{error.sd}, the standard deviations \code{error.sd * par["sd_rel"]} are used. This makes it possible to estimate a common factor of a given variance-covariance structure.
}

\value{
Returns log likelihood.
}

\seealso{
\code{\link{Logposterior}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}