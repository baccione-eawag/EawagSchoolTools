\name{Plot.chains}

\alias{Plot.chains}

\title{Plotting functions (EawagSchoolTools)}

\description{
Function to plot components of a Markov chain.
}

\usage{
Plot.chains(postsamp,ncol=NA,mar=NA,ylim=list(),titles=list(),
            xlab="chain index",ylab=list())
}

\arguments{
  \item{postsamp}{data frame containing Markov chain}
  \item{ncol}{optional number of columns in plotting grid}
  \item{mar}{optional margins}
  \item{ylim}{optional limits of the y axes (list)}
  \item{titles}{optional titles (list)}
  \item{xlab}{optional labels of the x axes}
  \item{ylab}{optional labels of the y axes (list)}
}

\value{
Plots of components of Markov chain.
}

\seealso{
\code{\link{Plot.res}}, \code{\link{Plot.vars}}, \code{\link{contourpdf}}, \code{\link{Plot.margs}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}