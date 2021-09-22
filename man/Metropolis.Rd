\name{Markovchain.Metropolis}

\alias{Markovchain.Metropolis}

\title{Metropolis Markov chain sample}

\description{
Calculates a Markov Chain of a probability distribution of a vector of continuous random variables, using the Metropolis algorithm with a normal jump distribution with given standard deviations and correlation structure.
}

\usage{
Markovchain.Metropolis(log.pdf,z.ini,prop.sd,prop.cor=0,
                       sampsize,thin=1,...)
}

\arguments{
  \item{log.pdf}{function \code{log.pdf(z,...)} that calculates the log of the probability density function of the vector of random variables that are to be sampled. \code{log.pdf} must be given a single vector \code{z} at which the density is to be evaluated. Additional arguments from the call to \code{markovchain.metropolis} will be passed on. If the probability density is zero, NA must be returned.}
  \item{z.ini}{vector of values at which the chain is to be started}
  \item{prop.sd}{vector of standard deviations of the proposal distribution}
  \item{prop.cor}{correlation matrix of proposal distribution, or NA if all correlations are zero}
  \item{sampsize}{sample size (length of the chain)}
  \item{thin}{factor with which to thin storage of results (thin=n: only each nth result is returned; this saves memory)}
  \item{...}{additional parameters passed to \code{log.pdf}}
}

\details{
The log of the probability density of the distribution must be specified as a function \code{log.pdf(z,...)}, where \code{z} is the vector of values for which the density should be evaluated.
}

\value{
Returns a list containing:
  \item{z}{sample as a matrix, with sample points in rows}
  \item{log.pdf}{vector with log pdf values of the sample}
  \item{reject.freq}{rejection frequency of the jumps}
  \item{error}{error message (empty string if no error occurred)}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}