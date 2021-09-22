################################################################################
#                                                                              #
# Functions for R package "bayesBias" : Bayesian Inference with Consideration  #
# of Bias                                                                      #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# Function to calculate probability densities of univariate distributions:
# ------------------------------------------------------------------------

calcpdf <- function(x,distpar,log=FALSE)
{
  if ( distpar[1] == "Uniform" | distpar[1] == "uniform" )
  {
    # uniform distribution; parameters are min and max
    min <- as.numeric(distpar[2])
    max <- as.numeric(distpar[3])
    return(dunif(x,min=min,max=max,log=log))
  }
  if ( distpar[1] == "Normal" | distpar[1] == "normal" )
  {
    # normal distribution; parameters are mean and sd:
    mean <- as.numeric(distpar[2])
    sd   <- as.numeric(distpar[3])
    return(dnorm(x,mean=mean,sd=sd,log=log))
  }
  if ( distpar[1] == "NormalTrunc" | distpar[1] == "normaltrunc" )
  {
    # truncated normal distribution; parameters are mean, sd, min and max
    mean <- as.numeric(distpar[2])
    sd   <- as.numeric(distpar[3])
    min  <- as.numeric(distpar[4])
    max  <- as.numeric(distpar[5])
    fact <- 1/(pnorm(q=max,mean=mean,sd=sd)-pnorm(q=min,mean=mean,sd=sd))
    if ( !log )
    {
      return(ifelse(x<min|x>max,0,fact*dnorm(x,mean=mean,sd=sd)))
    }
    else
    {
      return(ifelse(x<min|x>max,NA,
                    log(fact)+dnorm(x,mean=mean,sd=sd,log=TRUE)))
    }
  }
  if ( distpar[1] == "Lognormal" | distpar[1] == "lognormal")
  {
    # lognormal distribution; parameters are mean and sd;
    # R parameters are mean and sd of the log of the random variable
    mean    <- as.numeric(distpar[2])
    sd      <- as.numeric(distpar[3])
    sdlog   <- sqrt(log(1+sd^2/mean^2))
    meanlog <- log(mean) - 0.5*sdlog^2
    return(dlnorm(x,meanlog=meanlog,sdlog=sdlog,log=log))
  }
  if ( distpar[1] == "LognormalTrunc" | distpar[1] == "lognormaltrunc" )
  {
    # truncated lognormal distribution; parameters are mean, sd, min and max;
    # R parameters are mean and sd of the log of the random variable
    mean    <- as.numeric(distpar[2])
    sd      <- as.numeric(distpar[3])
    sdlog   <- sqrt(log(1+sd^2/mean^2))
    meanlog <- log(mean) - 0.5*sdlog^2
    min     <- as.numeric(distpar[4])
    max     <- as.numeric(distpar[5])
    fact <- 1/(plnorm(q=max,meanlog=meanlog,sdlog=sdlog)-plnorm(q=min,meanlog=meanlog,sdlog=sdlog))
    if ( !log )
    {
      return(ifelse(x<min|x>max,0,fact*dlnorm(x,meanlog=meanlog,sdlog=sdlog)))
    }
    else
    {
      return(ifelse(x<min|x>max,NA,
                    log(fact)+dlnorm(x,meanlog=meanlog,sdlog=sdlog,log=TRUE)))
    }
  }
  if ( distpar[1] == "Inv" | distpar[1] == "inv" )
  {
    # inverse distribution; parameters are min and max:
    min <- as.numeric(distpar[2])
    max <- as.numeric(distpar[3])
    if ( !log )
    {
      return(ifelse(x<min|x>max,0,1/(log(max/min)*x)))   
    }
    else
    {
      return(ifelse(x<min|x>max,NA,-log(log(max/min)) - log(x)))   
    }
  }
  if ( distpar[1] == "Exponential" | distpar[1] == "exponential" )
  {
    # exponential distribution; parameter is mean:
    mean <- as.numeric(distpar[2])
    if ( !log )
    {
      return(ifelse(x<0,0,1/mean*exp(-x/mean)))   
    }
    else
    {
      return(ifelse(x<0,NA,-log(mean)-x/mean))   
    }
  }
  stop(paste("Distribution",dist,"not yet implemented"))
}


################################################################################