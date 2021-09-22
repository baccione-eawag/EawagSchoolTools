################################################################################
#                                                                              #
# Functions for R package "bayesBias" : Bayesian Inference with Consideration  #
# of Bias                                                                      #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# Calculate the logarithm of the probability density function of a multivariate
# normal or lognormal distribution or of a product of independent marginals.

calcpdf_mv <- function(z,dist="normal",mean=0,sd=1,cor=0,
                          cor.inv=NA,log=TRUE,distdef=NA,file=NA)
{
  # -----------------------------------------------------------------------
  # This function calculates the logarithm of the probability density 
  # function of a multivariate normal or lognormal distribution or of 
  # independent 1d distributions.
  #
  # Arguments:
  # z:          vector, matrix or data frame at which the logarithm of the 
  #             probability density function has to be evaluated
  # dist:       distribution type: "Normal", "Lognormal" or "Indep".
  # mean:       vector of means
  # sd:         vector of standard deviations
  # cor:        correlation matrix of the distribution
  # cor.inv:    inverse of correlation matrix of the distribution
  #             (alternative input to cor, saves computation time for
  #             repeated calls)
  # log:        if TRUE returns log of pdf, otherwise pdf
  # distdef:    distribution definition for independent 1d distributions
  #
  # Return Value:
  # probability density (or log of) for all sample points in x
  #
  #                                        Peter Reichert    Jan.  01, 2004
  #                                        last modification Oct.  29, 2011
  # -----------------------------------------------------------------------
  
  # consistency checks and initializations:
  mean <- as.vector(mean)
  sd <- as.vector(sd)
  if ( length(sd) != length(mean) )
  {
    stop("calcpdf_mv: illegal dimension of standard deviations")
  }
  if ( is.vector(z) )
  {
    len <- length(z)
    names <- names(z)
    z <- as.matrix(z)
    dim(z) <- c(1,len)
    if ( length(names) == len ) colnames(z) <- names
  }
  numpar <- ncol(z)
  R <- diag(rep(1,numpar))
  if ( is.matrix(cor) ) R <- cor
  if ( nrow(R) != numpar || ncol(R) != numpar )
  {
    stop("calcpdf_mv: illegal dimension of correlation matrix")
  }
  
  # calculate logarithm of probability density function:
  sampsize <- nrow(z)
  logpdf <- numeric(sampsize)
  if ( dist == "normal" | dist == "Normal" )
  {
    # multivariate normal distribution:
    n <- length(sd)
    sigma <- diag(sd,nrow=n,ncol=n) %*% R %*% diag(sd,nrow=n,ncol=n)
    if ( is.matrix(cor.inv) )
    {
      R.inv <- cor.inv
    }
    else
    {
      R.inv <- solve(R)
    }
    det.R.inv <- det(R.inv)
    if ( det.R.inv > 0 )
    {
      for ( i in 1:sampsize )
      {
        v <- as.matrix(z[i,]-mean,nrow=numpar,ncol=1)/
          as.matrix(sd,nrow=numpar,ncol=1)
        logpdf[i] <- -numpar/2*log(2*pi) + 0.5*log(det.R.inv) - log(prod(sd)) -
          0.5 * t(v) %*% R.inv %*% v
      }
    }
    else
    {
      logpdf <- rep(NA,sampsize)
    }
  }
  else
  {
    if ( dist == "lognormal" | dist == "Lognormal" )
    {
      # multivariate lognormal distribution:
      sdlog    <- sqrt(log(1+sd*sd/(mean*mean)))
      meanlog  <- log(mean) - sdlog*sdlog/2
      if ( numpar > 1 )
      {
        n <- length(sdlog)
        ln.sigma <- log( 1 + diag(sqrt(exp(sdlog*sdlog)-1),nrow=n,ncol=n) %*% 
                           R %*% diag(sqrt(exp(sdlog*sdlog)-1),nrow=n,ncol=n) )
      }
      else
      {
        ln.sigma <- as.matrix( log( 1 + sqrt(exp(sdlog*sdlog)-1)^2 ) )
      }
      ln.sigma.inv = solve(ln.sigma)
      det.ln.sigma = det(ln.sigma)
      if ( det.ln.sigma > 0 )
      {
        for ( i in 1:sampsize )
        {
          if ( min(z[i,]) <= 0 )
          {
            logpdf[i] <- NA 
          }
          else
          {
            v <- log(z[i,])-meanlog
            v <- matrix(v,nrow=numpar,ncol=1) # ensure that v is a column vector
            logpdf[i] <- -numpar/2*log(2*pi) - 0.5*log(det.ln.sigma) - 
              log(prod(z[i,])) - 0.5 * t(v) %*% ln.sigma.inv %*% v
          }
        }
      }
      else
      {
        logpdf <- rep(NA,sampsize)
      }
    }
    else
    {
      if ( dist == "indep" | dist == "Indep" ) 
      {
        for ( i in 1:sampsize )
        {
          logpdf[i] <- 0
          for ( j in 1:ncol(z) )
          {
            ind <- j
            if ( length(colnames(z)) == ncol(z) )
            {
              ind <- match(colnames(z)[j],names(distdef))
              if ( is.na(ind) )
              {
                stop(paste("error in calcpdf_mv:",
                           "variable",colnames(z)[j],"not found"))
              }
            }
            logpdf[i] <- logpdf[i] + calcpdf(z[i,j],
                                             distdef[[ind]],
                                             log=TRUE)
          }
        }
      }
      else
      {
        stop(paste("calcpdf_mv: unknown distribution type:",dist))
      }
    }
  }
  
  # write results:
  if ( !is.na(file) )
  {
    write.table(data.frame(z,logpdf=logpdf),file=file,
                col.names=TRUE,row.names=FALSE,sep="\t")
  }
  
  # return result:
  if ( log )
  {
    return(logpdf)
  }
  else
  {
    return(exp(logpdf))
  }
}


################################################################################