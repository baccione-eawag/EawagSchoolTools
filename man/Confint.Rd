\name{Confint}

\alias{Confint}

\title{Confidence intervals (eawagSummerSchoolTools)}

\description{
Calculate confidence intervals based on the Student t-distribution or on the normal distribution.
}

\usage{
Confint(est,sd,df,level)
}

\arguments{
  \item{est}{vector of point estimates}
  \item{sd}{vector of estimated standard deviations}
  \item{df}{degrees of freedom}
  \item{level}{probability level of confidence intervals (default 0.95)}
}

\details{
Confidence intervals are calculated based on the normal distribution if \code{df=NA}.
}

\value{
Returns a table of lower and upper bounds of confidence intervals.
}

\seealso{
\code{\link{Gnls}}, \code{\link{Gnls.diag}}, \code{\link{Gnls.test}}, \code{\link{Gnls.predict}}, \code{\link{Confint}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}