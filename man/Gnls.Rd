\name{Gnls}

\alias{Gnls}

\title{Generalized nonlinear least squares parameter estimates}

\description{
Calculate the generalized least squares parameter estimates of a nonlinear model.
}

\usage{
Gnls(model.name,y.obs,var.inv,par,x=list(),...)
}

\arguments{
  \item{model.name}{name of the function representing the model (note that this model requires the parameters to be specified explicitly in the function headers)}
  \item{y.obs}{observed data corresponding to model output}
  \item{var.inv}{inverse variance-covariance matrix of the multivariate normal distribution that characterizes the error model}
  \item{par}{named parameter vector with initial values passed to the model as separate arguments}
  \item{x}{list of named inputs passed to the model}
  \item{...}{additional arguments passed to \code{nls}}
}

\value{
Returns a list containing:
  \item{model.name}{name of the function representing the model}
  \item{par}{named parameter vector with parameter estimates}
  \item{y.obs}{observed data corresponding to model output}
  \item{y.det}{deterministic model results corresponding to the estimated parameters}
  \item{resid}{residuals}
  \item{var.inv}{inverse variance-covariance matrix of the multivariate normal distribution that characterizes the error model}
  \item{x}{input object passed as first argument to the model}
  \item{res.nls}{results from application of the \code{nls} function}
}

\seealso{
\code{\link{Gnls}}, \code{\link{Gnls.diag}}, \code{\link{Gnls.test}}, \code{\link{Gnls.predict}}, \code{\link{Confint}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}