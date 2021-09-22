\name{Plot.margs}

\alias{Plot.margs}

\title{Plotting functions (eawagSummerSchoolTools)}

\description{
Function to plot marginal distributions of a sample. If definitions of prior marginals are provided, these can be plotted in addition to the sample marginals.
}

\usage{
Plot.margs(postsamp,pridist=list(),val=NA,ncol=NA,nrow=NA,mar=NA,
           xlim=list(),ymax=list(),titles=list(),xlab=list(),ylab=list(),
           adjust=1,lty=NA,col.post=NA,col.pri=NA)
}

\arguments{
  \item{postsamp}{data frame containing Markov chain}
  \item{pridist}{prior distributions (list)}
  \item{val}{optional values}
  \item{ncol}{optional number of columns in plotting grid}
  \item{nrow}{optional number of rows in plotting grid}
  \item{mar}{optional margins}
  \item{xlim}{optional limits of the x axes (list)}
  \item{ymax}{optional maximum y values (list)}
  \item{titles}{optional titles (list)}
  \item{xlab}{optional labels of x axes (list)}
  \item{ylab}{optional labels of y axes (list)}
  \item{adjust}{optional adjustment}
  \item{lty}{optional line types}
  \item{col.post}{optional color for posteriors}
  \item{col.pri}{optional color for priors}
}

\value{
Plots of marginal distributions.
}

\seealso{
\code{\link{Plot.res}}, \code{\link{Plot.vars}}, \code{\link{contourpdf}}, \code{\link{Plot.chains}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}