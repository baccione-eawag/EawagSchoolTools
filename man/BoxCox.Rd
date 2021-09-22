\name{BoxCox}

\alias{BoxCox}
\alias{BoxCox.deriv}
\alias{BoxCox.inv}

\title{Box-Cox transformation functions}

\description{
Functions to calculate the Box-Cox transformation, its derivative, and inverse transformation (back-transformation).
}

\usage{
BoxCox(      data,lambda1=1,lambda2=1)
BoxCox.deriv(data,lambda1=1,lambda2=1)
BoxCox.inv(  data,lambda1=1,lambda2=1)
}

\arguments{
  \item{data}{
    data to be transformed (\code{boxcox}, \code{boxcox.deriv}) or back-transformed (\code{boxcox.inv})
  }
  \item{lambda1}{
    value of parameter lambda1 (exponent) of Box-Cox transformation
  }
  \item{lambda2}{
    value of parameter lambda2 (offset) of Box-Cox transformation
  }
}

\value{
Transformed data (\code{boxcox}), derivative of Box-Cox transformation (\code{boxcox.deriv}), or back-transformed data (\code{boxcox.inv}).
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}