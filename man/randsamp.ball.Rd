\name{randsamp.ball}

\alias{randsamp.ball}

\title{Random sample from an n-dimensional ball}

\description{
Generate a random sample from a uniform distribution in an n-dimensional ball of given radius.
}

\usage{
randsamp.ball(sampsize=1,dim=1,radius=1)
}

\arguments{
  \item{sampsize}{sample size}
  \item{dim}{dimension of ball}
  \item{radius}{radius of ball}
}

\details{
For details and references, see introduction of: Harman, R. and Lacko, V. (2010), \emph{On decompositional algorithms for uniform sampling from n-spheres and n-balls}, Journal of Multivariate Analysis 101, 2297-2304.
}

\value{
Sampling output.
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}