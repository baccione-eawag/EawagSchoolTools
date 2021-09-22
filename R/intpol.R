################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


intpol.irregular <- function(x,
                             x.data,
                             y.data,               
                             method=c("neighbor","invdist"),
                             scales.x=NA,
                             fact=1)
{
  #   interpolation on an irregular grid
  #
  #   x:        vector of values of independent variable at which to interpolate;
  #             xout can also be a matrix, interpolation is then for each row
  #   x.data:   matrix of independent independent variables; each row of x.data
  #             belongs to one input specification (x.data can be a vector if there
  #             is only one input dimension)
  #   y.data:   vector of values of the dependent variable; one value for each
  #             row of x.data
  #   method:   "neighbor" does nearest neighbor interpolation (piecewise 
  #             constant function of value at nearest data point)
  #             "invdist2" weights the provided dependent values according to 
  #             the inverse of the distance to xout squared
  #   scales.x: values of x and of x.data are scaled by these scales before 
  #             interpolation; default is all scales equal to 1
  #
  #   --------------------------------------
  
  
  # call of each row if x is a matrix:
  
  if ( is.vector(x.data) )
  {
    x.data <- matrix(x.data,ncol=1)
    x      <- matrix(x,ncol=1)
  }
  if ( is.matrix(x) )
  {
    y <- rep(NA,nrow(x))
    for ( i in 1:nrow(x) )
    {
      y[i] <- intpol.irregular(as.vector(x[i,]),x.data,y.data,method,scales.x,fact)
    }
    return(y)
  }
  
  # check input:
  
  n.data <- length(y.data)
  dim.x  <- length(x)
  y <- NA
  if ( nrow(x.data)!=n.data | ncol(x.data)!=dim.x )
  {
    print("incorrect dimensions in interpol.irregular")
    return(y)
  }
  if ( is.na(scales.x[1]) ) scales.x <- rep(1,dim.x)
  
  # calculate distances to provided points:
  
  x.mat <- matrix(rep(x,n.data),nrow=n.data,byrow=TRUE)
  dist <- sqrt(apply(((x.mat-x.data)%*%diag(1/scales.x))^2,1,sum))
  
  ind.min <- which.min(dist)
  if ( method[1] == "neighbor" )
  {
    y <- y.data[ind.min]
  }
  else
  {
    if ( method[1] == "invdist" )
    {
      if ( dist[ind.min] == 0)
      {
        y <- y.data[ind.min]
      }
      else
      {
        w <- 1/dist^(dim.x*fact)
        w <- w/sum(w)
        y <- sum(w*y.data)
      }
    }
    else
    {
      print(paste("unknown method:",method[1]))
    }
  }
  
  return(y)
}


################################################################################