\name{Plot.vars}

\alias{Plot.vars}

\title{Plotting functions (eawagSummerSchoolTools)}

\description{
Function to plot several sets of results.
}

\usage{
Plot.vars(vars,L,ncol=NA,mar=NA,ylim=list(),markers=F,headers=list(),
          xlab="",ylab="",pos="topright")
}

\arguments{
  \item{vars}{matrix or data frame with variables, and variable information encoded in row names}
  \item{L}{observation layout}
  \item{ncol}{optional number of columns of sub-panels of the plot}
  \item{mar}{optional specification of margins in the form c(bottom,left,top,right)}
  \item{ylim}{optional named (by variable name) list of limits of the y axes}
  \item{markers}{if TRUE plot markers instead of lines}
  \item{headers}{optional named (by variable name) list of headers of the plots}
  \item{xlab}{optional label of the x axis}
  \item{ylab}{optional label of the y axis}
  \item{pos}{position of legend (only if more than one variable)}
}

\details{
Results should be stored in a data frame or matrix, with the corresponding variable information encoded in the row names (for example C_M_1 for variable C_M at time 1).
}

\value{
Plot of all variables as a function of the independent variable.
}

\seealso{
\code{\link{Plot.res}}, \code{\link{contourpdf}}, \code{\link{Plot.chains}}, \code{\link{Plot.margs}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}