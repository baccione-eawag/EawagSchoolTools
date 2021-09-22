\name{Hessian}

\alias{Hessian}

\title{Hessian}

\description{
Calculates the Hessian of a real-valued function.
}

\usage{
Hessian(fn,par,par.inc=NA,...)
}

\arguments{
  \item{fn}{function}
  \item{par}{parameters of the function}
  \item{par.inc}{optional increments}
  \item{...}{additional arguments passed to \code{fn}}
}

\value{
Returns a list containing:
  \item{h}{hessian matrix}
  \item{minimum}{TRUE or FALSE}
  \item{maximum}{TRUE or FALSE}
  \item{num.eval}{number of evaluations}
  \item{trace}{trace (matrix)}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}