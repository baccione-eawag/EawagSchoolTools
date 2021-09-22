################################################################################
#                                                                              #
# Functions for R package "bayesBias" : Bayesian Inference with Consideration  #
# of Bias                                                                      #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# Generate a random sample from a multivariate normal or lognormal distribution 
# or from a product of independent 1d marginals.

randsamp <- function(sampsize=1,dist="Normal",mean=0,sd=1,cor=NA,
                     distdef=NA,file=NA)
{
  # -----------------------------------------------------------------------
  # This function generates a random sample from a multivariate
  # normal or lognormal distribution.
  # This function is a simplified version of the program "randsamp"
  # available as part of the package UNCSIM at http://www.uncsim.eawag.ch
  #
  # Arguments:
  # sampsize:   sample size
  # mean:       vector of means
  # sd:         vector of standard deviations
  # cor:        correlation matrix of the distribution
  # dist:       distribution type: "Normal", "Lognormal" or "Indep".
  # distdef:    definition of 1d marginals for "Indep" instead of (mean,sd,cor)
  #
  # Return Value:
  # List of:
  # mean:       vector of means
  # sd:         vector of standard deviations
  # corr:       correlation matrix of the distribution
  # sampsize:   sample size
  # sample:     data frame of parameter samples (each row corresponds 
  #             to a draw)
  # logpdf:     vector of values of log probability density at the sample
  #             points
  #
  #                                        Peter Reichert    Dec.  29, 2003
  #                                        last modification Oct.  29, 2011
  # -----------------------------------------------------------------------
  
  R <- 0
  if ( dist == "normal"    | dist == "Normal" | 
         dist == "lognormal" | dist == "Lognormal")
  {
    # consistency checks and initializations:
    
    numpar <- length(mean)
    if ( length(sd) != numpar )
    {
      stop("randsamp: mean and sd do not have the same length")
    }
    R <- diag(rep(1,numpar))
    if ( is.matrix(cor) ) R <- cor
    if ( nrow(R) != numpar || ncol(R) != numpar )
    {
      stop("randsamp: illegal dimension of correlation matrix")
    }
    
    # calculate sample from multivariate uniform distribution:
    samp <- runif(sampsize*numpar)
    dim(samp) <- c(sampsize,numpar)
    
    # transform sample to multivariate normal or lognormal and calculate
    # logarithm of probability density function:
    logpdf <- numeric(sampsize)
    if ( dist == "normal" | dist == "Normal" )
    {
      # calculate transformation matrix and transform the sample:
      sigma <- diag(sd) %*% R %*% diag(sd)
      sigma.inv = solve(sigma)
      det.sigma = det(sigma)
      A <- t(chol(sigma))
      for ( i in 1:sampsize )
      {
        samp[i,]  <- A %*% qnorm(samp[i,]) + mean
        v <- samp[i,]-mean
        v <- matrix(v,nrow=numpar,ncol=1) # ensure that v is a column vector
        logpdf[i] <- -numpar/2*log(2*pi) - 0.5*log(det.sigma) -
          0.5 * t(v) %*% sigma.inv %*% (v)
      }
    }
    else # dist == "lognormal" | dist == "Lognormal"
    {
      # parameters of the log of the variable, calculate transformation 
      # matrix and transform the sample:
      sdlog    <- sqrt(log(1+sd*sd/(mean*mean)))
      meanlog  <- log(mean) - sdlog*sdlog/2
      if ( numpar > 1 )
      {
        ln.sigma <- log( 1 + diag(sqrt(exp(sdlog*sdlog)-1)) %*%
                           R %*% diag(sqrt(exp(sdlog*sdlog)-1)) )
      }
      else
      {
        ln.sigma <- as.matrix( log( 1 + sqrt(exp(sdlog*sdlog)-1)^2 ) )
      }
      ln.sigma.inv = solve(ln.sigma)
      det.ln.sigma = det(ln.sigma)
      ln.A <- t(chol(ln.sigma))
      for ( i in 1:sampsize )
      {
        log.samp.i <- ln.A %*% qnorm(samp[i,]) + meanlog
        samp[i,]  <- exp(log.samp.i)
        v <- log.samp.i-meanlog
        v <- matrix(v,nrow=numpar,ncol=1) # ensure that v is a column vector
        logpdf[i] <- -numpar/2*log(2*pi) - 0.5*log(det.ln.sigma) - 
          log(prod(samp[i,])) - 0.5 * t(v) %*% ln.sigma.inv %*% v
      }
    }
    
    # collect results:
    colnames(samp) <- names(mean)
    samp <- data.frame(samp)
    res <- list(
      mean       = mean,
      sd         = sd,
      cor        = R,
      sampsize   = sampsize,
      sample     = samp,
      logsamppdf = logpdf
    )
  }
  else   # dist != normal,Normal,lognormal,Lognormal      
  {  
    if ( dist == "indep" | dist == "Indep" )
    {
      logpdf <- 0
      samp <- matrix(NA,nrow=sampsize,ncol=length(distdef))
      colnames(samp) <- names(distdef)
      samp <- data.frame(samp)
      for ( j in 1:ncol(samp) )
      {
        distpar <- distdef[[j]]
        dist.found <- F
        if ( !dist.found )
        {
          if ( distpar[1] == "Uniform" | 
                 distpar[1] == "uniform" )
          {
            # uniform distribution; parameters are min and max
            min <- as.numeric(distpar[2])
            max <- as.numeric(distpar[3])
            samp[,j] <- runif(sampsize,min=min,max=max)
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "Normal" | 
                 distpar[1] == "normal" )
          {
            # normal distribution; parameters are mean and sd:
            mean <- as.numeric(distpar[2])
            sd   <- as.numeric(distpar[3])
            samp[,j] <- rnorm(sampsize,mean=mean,sd=sd)
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "NormalTrunc" | 
                 distpar[1] == "normaltrunc" )
          {
            # truncated normal distribution; parameters are mean, sd, min and max
            # of untruncated normal distribution
            mean <- as.numeric(distpar[2])
            sd   <- as.numeric(distpar[3])
            min  <- as.numeric(distpar[4])
            max  <- as.numeric(distpar[5])
            cdf.min <- pnorm(min,mean=mean,sd=sd)
            cdf.max <- pnorm(max,mean=mean,sd=sd)
            samp[,j] <- runif(sampsize,min=cdf.min,max=cdf.max)
            samp[,j] <- qnorm(samp[,j],mean=mean,sd=sd)
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "Lognormal" | 
                 distpar[1] == "lognormal" )
          {
            # lognormal distribution; parameters are mean and sd:
            mean    <- as.numeric(distpar[2])
            sd      <- as.numeric(distpar[3])
            sdlog   <- sqrt(log(1+sd^2/mean^2))
            meanlog <- log(mean) - 0.5*sdlog^2
            samp[,j] <- rlnorm(sampsize,meanlog=meanlog,sdlog=sdlog)
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "LognormalTrunc" | 
                 distpar[1] == "lognormaltrunc" )
          {
            # truncated lognormal distribution; parameters are mean, sd, min and max
            # of untruncated lognormal distribution
            mean    <- as.numeric(distpar[2])
            sd      <- as.numeric(distpar[3])
            sdlog   <- sqrt(log(1+sd^2/mean^2))
            meanlog <- log(mean) - 0.5*sdlog^2
            min     <- as.numeric(distpar[4])
            max     <- as.numeric(distpar[5])
            cdf.min <- plnorm(min,meanlog=meanlog,sdlog=sdlog)
            cdf.max <- plnorm(max,meanlog=meanlog,sdlog=sdlog)
            samp[,j] <- runif(sampsize,min=cdf.min,max=cdf.max)
            samp[,j] <- qlnorm(samp[,j],meanlog=meanlog,sdlog=sdlog)
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "Inv" | 
                 distpar[1] == "inv" )
          {
            # inverse distribution (f(x) prop. 1/x); 
            # parameters are min and max:
            min     <- as.numeric(distpar[2])
            max     <- as.numeric(distpar[3])
            log.min <- log(min)
            log.max <- log(max)
            samp[,j] <- runif(sampsize,min=0,max=1)
            samp[,j] <- exp(samp[,j]*(log.max-log.min)+log.min)
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "Exponential" | 
                 distpar[1] == "exponential" )
          {
            # exponential distribution; parameter is mean:
            mean <- as.numeric(distpar[2])
            samp[,j] <- runif(sampsize,min=0,max=1)
            samp[,j] <- -mean * log(1-samp[,j])
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "Delta" | 
                 distpar[1] == "delta" )
          {
            # delta distribution; parameter is mean:
            mean <- as.numeric(distpar[2])
            samp[,j] <- rep(sampsize,mean)
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "Discrete" | 
                 distpar[1] == "discrete" )
          {
            # discrete distribution; parameters are probabilities:
            n <- length(distpar[-1])/2
            probs <- as.numeric(distpar[1+n+(1:n)])
            probs <- probs/sum(probs)
            probs.lowerbounds <- rep(NA,n)
            probs.lowerbounds[1] <- 0
            if ( n > 1 )
            {
              for ( i in 2:n ) probs.lowerbounds[i] <- 
                probs.lowerbounds[i-1] + probs[i-1]
            }
            samp.unif <- runif(sampsize,min=0,max=1)
            ind <- rep(NA,sampsize)
            for ( i in 1:sampsize )
            {
              ind[i] <- sum(samp.unif[i] > probs.lowerbounds)
            }
            samp[,j] <- distpar[1+ind]
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          stop(paste("Distribution",distpar[1],"not yet implemented"))
        }
      }
    }
    else
    {
      stop(paste("randsamp: unknown distribution type:",dist))
    }
  }
  
  # collect results:
  res <- list(
    mean       = mean,
    sd         = sd,
    cor        = R,
    sampsize   = sampsize,
    sample     = samp,
    logsamppdf = logpdf
  )
  
  # write results:
  if ( !is.na(file) )
  {
    write.table(data.frame(samp,logsamppdf=logpdf),file=file,
                col.names=TRUE,row.names=FALSE,sep="\t")
  }
  
  # return results:
  return(res)
}


################################################################################