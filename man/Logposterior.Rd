\name{Logposterior}

\alias{Logposterior}

\title{Log posterior}

\description{
Function to calculate the log posterior probability density for a model. It is based on a likelihood definition provided by the user, and a multivariate normal or lognormal distribution at the prior.

This log posterior implementation serves as a template for implementing more specific posteriors.
}

\usage{
Logposterior(par,model,L,y,prior.dist="lognormal",
             prior.mean=1,prior.sd=1,prior.cor=NA,prior.def=NA,
             loglikeli=Loglikeli,...)
}

\arguments{
  \item{par}{parameters}
  \item{model}{model function}
  \item{L}{observation layout}
  \item{y}{data}
  \item{prior.dist}{prior distribution type (default "lognormal")}
  \item{prior.mean}{prior distribution mean}
  \item{prior.sd}{prior distribution standard deviation}
  \item{prior.cor}{prior distribution correlation matrix}
  \item{prior.def}{prior distribution definition (if using independent one-dimensional distributions)}
  \item{loglikeli}{log likelihood function}
  \item{...}{additional arguments passed to \code{loglikeli}}
}

\value{
Returns log posterior.
}

\seealso{
\code{\link{Loglikeli}}.
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}