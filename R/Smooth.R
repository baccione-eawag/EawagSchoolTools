################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# Smooth:
# -------
#
# Function for smoothing data points and estimating the derivative of the 
# smoothed curve by local quadratic or optionally local linear regression. 
# Local regression is implemented by using a Gaussian distribution of 
# of weights centered a the point at which the smoothed curve has to be
# evaluated.
# The smoothing parameter is the standard deviation of the Gaussian weights.
#
# Peter Reichert 05.09.2008 , last modification 01.05.2009


# Call:      Smooth(x,y,sigma,newx=NA,fac.extrap=1,method="quadratic")
# -----
#
# Input:
# ------
#
# x          vector of x-coordinates of data points to be smoothed
# y          vector of y-coordinates of data points to be smoothed
#            (x and y must be of the same length)
# sigma      standard deviation of Gaussian distribution used as weights for
#            local quadratic or local linear regression
# newx       optional vector of x-coordinates at which smoothed results and
#            derivatives are to be calculated (if not specified, results
#            are provided at the same locations as there is data available)
# fac.extrap calculate smoothed value only if data is available within
#            fac.extrap*sigma (default is 1)
# method     "quadratic" (default) indicates local quadratic regression, 
#            otherwise local linear regression is used
#
# Output:
# -------
#
# Data frame of
# x          x-coordinates at which smoothed results and derivatives are 
#            available
# y          smoothed results at the locations x
# ydot       derivatives of smoothed results at the locations x

Smooth <- function(x,y,sigma,newx=NA,fac.extrap=1,method="quadratic")
{
  # calculate x and y vectors for available data:
  ind <- !is.na(y)
  x.data <- x[ind]
  y.data <- y[ind]
  
  # check consistency of input:
  if ( length(x.data) < 1 )
  {
    stop("*** error in Smooth: no data available ***")
  }
  if ( length(x.data) != length(y.data) ) 
  {
    stop("*** error in Smooth: length of x and y different ***")
  }
  if ( ! sigma > 0 ) 
  {
    stop("*** error in Smooth: sigma is not positive ***")
  }
  
  # select x values for output:
  if ( is.na(newx[1]) ) newx <- x
  
  # set up ysmooth and ydot vectors:
  n <- length(newx)
  ysmooth <- rep(NA,n)
  ydot    <- rep(NA,n)
  
  # calclate smoothed values and derivatives:
  for ( i in 1:n )
  {
    # get indices of data within +/- 2*sigma:
    ind.extrap <- x.data >= newx[i]-fac.extrap*sigma & 
      x.data <= newx[i]+fac.extrap*sigma
    num.extrap <- sum(ifelse(ind.extrap,1,0))
    
    # calculate smoothed value only if data is available within +/- 2*sigma:
    if ( num.extrap > 0 ) 
    {
      # still use data within a 5 times larger interval
      # to calculate the smoothed value:
      fac.use <- 4*max(1,fac.extrap)
      ind.use <- x.data >= newx[i]-fac.use*sigma & 
        x.data <= newx[i]+fac.use*sigma
      x1  <- x.data[ind.use]-newx[i]
      x2  <- (x.data[ind.use]-newx[i])^2
      num.use <- sum(ifelse(ind.use,1,0))
      if ( num.use == 1 )  # use value
      {
        ysmooth[i] <- y.data[ind.use][1]
        ydot[i]    <- 0
      }
      else
      {
        if ( num.use == 2 )  # use weighted mean
        {
          weights <- dnorm(x.data[ind.use],mean=newx[i],sd=sigma)
          weights <- weights/sum(weights)
          ysmooth[i] <- weights[1]*y.data[ind.use][1] + 
            weights[2]*y.data[ind.use][2]
          if ( x.data[ind.use][2] != x.data[ind.use][1] )
          {
            ydot[i]    <- (y.data[ind.use][2]-y.data[ind.use][1])/
              (x.data[ind.use][2]-x.data[ind.use][1])
          }
        }
        else
        {
          if ( method != "quadratic" | num.use == 3 ) # use local linear
          {                                           # regression
            res.lm     <- lm(y.data[ind.use] ~ x1,
                             weights=dnorm(x.data[ind.use],
                                           mean=newx[i],sd=sigma))
            ysmooth[i] <- coef(res.lm)[1]
            ydot[i]    <- coef(res.lm)[2]
          }
          else  # use local quadratic regression
          {
            res.lm     <- lm(y.data[ind.use] ~ x1 + x2,
                             weights=dnorm(x.data[ind.use],
                                           mean=newx[i],sd=sigma))
            ysmooth[i] <- coef(res.lm)[1]
            ydot[i]    <- coef(res.lm)[2]
          }
        }
      }
    }
  }
  
  # return data frame:
  return(data.frame(x=newx,y=ysmooth,ydot=ydot))
}


################################################################################


# Smooth.fun:
# -----------
#
# Function for smoothing piecewise linear functions.
# The smoothing parameter is the standard deviation of the Gaussian weights.
#
# Peter Reichert 16.12.2008 , last modification 01.05.2009


# Call:      Smooth.fun(data,z,sigma,newz=NA,newx=NA,fac.extrap=1,
#                       method="quadratic")
# -----
#
# Input:
# ------
#
# data       list of matrices specifying piecewise linear functions:
#            for each function, the independent variable x must be provided 
#            in the first column, the dependent variable y in the second
# z          vector of z-coordinates corresponding to the functions
# sigma      standard deviation of Gaussian distribution used as weights for
#            local quadratic regression in z-coordinates (see Smooth)
# newz       optional vector of z-coordinates at which smoothed functions
#            are to be calculated (if not specified, results are provided
#            at the same locations as there is data available)
# newx       optional vector of x-coordinates  at which smoothed function 
#            values are to be calculated
# fac.extrap calculate smoothed value only if data is available within
#            fac.extrap*sigma (default is 1)
# method     "quadratic" (default) indicates local quadratic regression, 
#            otherwise local linear regression is used
#
# Output:
# -------
#
# List of data frames with columns
# x          x-coordinates at which smoothed results are available
# y          smoothed results at the locations x

Smooth.fun <- function(data,z,sigma,newz=NA,newx=NA,fac.extrap=1,
                       method="quadratic")
{
  # check consistency of input:
  n <- length(data)
  if ( length(z) != n )
  { 
    stop("error in Smooth_fun: not same number of locations as functions")
  }
  
  # select z and x values for output:
  if ( is.na(newz[1]) ) newz <- z
  if ( is.na(newx[1]) )
  {
    x <- numeric(0)
    for ( j in 1:n ) x <- c(x,data[[j]][,1])
    newx <- sort(unique(x))
  }
  
  # interpolate input functions to selected x values:    
  y <- matrix(nrow=length(newx),ncol=n)
  for ( j in 1:n )
  {
    y[,j] <- approx(x=data[[j]][,1],y=data[[j]][,2],xout=newx)$y
  }
  
  # set up result data structure:
  res <- list()
  for ( j in 1:length(newz) )
  {
    res[[j]] <- data.frame(x=newx,y=rep(NA,length(newx)))
  }
  names(res) <- newz
  
  # calculate results:
  for ( i in 1:length(newx) )
  {
    res.smooth <- Smooth(x=z,y=y[i,],sigma=sigma,newx=newz,
                         fac.extrap=fac.extrap,method=method)$y
    for ( j in 1:length(newz) )
    {
      res[[j]]$y[i] <- res.smooth[j] 
    }     
  }
  
  # return results:
  return(res)
}


################################################################################