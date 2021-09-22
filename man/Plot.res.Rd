\name{Plot.res}

\alias{Plot.res}

\title{Plotting functions (eawagSummerSchoolTools)}

\description{
Function to plot the results of a simulation.
}

\usage{
Plot.res(res,L,xlim=NA,ylim=NA,markers=F,header="",xlab="",ylab="",
         pos="topright")
}

\arguments{
  \item{res}{named vector of results}
  \item{L}{observation layout}
  \item{xlim}{optional limits of the x axis}
  \item{ylim}{optional limits of the y axis}
  \item{markers}{if TRUE plot markers instead of lines}
  \item{header}{optional header of the plot}
  \item{xlab}{optional label of the x axis}
  \item{ylab}{optional label of the y axis}
  \item{pos}{position of legend}
}

\details{
Results should be stored in a named numeric vector, with the corresponding variable information encoded in the names (for example C_M_1 for variable C_M at time 1).
}

\value{
Plot of all variables as a function of the independent variable.
}

\seealso{
\code{\link{Plot.vars}}, \code{\link{contourpdf}}, \code{\link{Plot.chains}}, \code{\link{Plot.margs}}
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}