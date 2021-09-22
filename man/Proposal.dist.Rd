\name{Proposal.dist}

\alias{Proposal.dist}

\title{Proposal distribution for a Metropolis Markov chain}

\description{
Function to calculate a proposal distribution for a Metropolis Markov chain. Provides the definitions of a multivariate normal jump distribution from a Markov chain sample.
}

\usage{
Proposal.dist(postsamp,fact.sd,fract.burnin=0,fact.cor=1,plot=F)
}

\arguments{
  \item{postsamp}{postsamp}
  \item{fact.sd}{fact.sd}
  \item{fract.burnin}{fract.burnin}
  \item{fact.cor}{fact.cor}
  \item{plot}{if TRUE returns a plot}
}

\value{
Returns a list containing:
  \item{sd}{standard deviation}
  \item{cor}{correlation}
}

\seealso{
\code{\link{Markovchain.Metropolis}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}