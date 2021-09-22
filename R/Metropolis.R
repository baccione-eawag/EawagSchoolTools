################################################################################
#                                                                              #
# Functions for R package "bayesBias" : Bayesian Inference with Consideration  #
# of Bias                                                                      #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# Calculate a Metropolis Markov Chain sample of a distribution

Markovchain.Metropolis <- function(log.pdf,z.ini,
                                   prop.sd,prop.cor=0,
                                   sampsize,thin=1,...)
{
  # -----------------------------------------------------------------------
  # This function calculates a Markov Chain of a probability distribution
  # of a vector of continuous random variables using the Metropolis
  # algorithm with a normal jump distribution with given standard
  # deviations and correlation structure.
  # The log of the probability density of the distribution must be 
  # specified as a function log.pdf(z,...) where z is the vector of values
  # for which the density has to be evaluated.
  #
  # Arguments:
  # log.pdf      function "log.pdf(z,...)" that calculates the log of the
  #              probability density function of the vector of random 
  #              variables that are to be sampled.
  #              log.pdf must be given a single vector z at which the
  #              density has to be evaluated. Additional arguments from
  #              the call to calc.markovchain.metropolis will be passed on.
  #              If the probability density is zero, NA must be returned.
  # z.ini        vector of values at which the chain is to be started.
  # prop.sd      vector of standard deviations of the jump distribution.
  # prop.cor     correlation matrix of jump distribution or NA if all
  #              correlations are zero.
  # sampsize     sample size (length of the chain)
  # thin         factor with which to thin storage of results (thin=n: 
  #              only each nth result is returned; this saves memory)
  # ...          additional parameters passed to "log.pdf"
  #
  # Return Value:
  # List with the following elements:
  # z            sample as a matrix with sample points in its rows.
  # log.pdf      vector with log pdf values of the sample.
  # reject.freq  rejection frequency of the jumps.
  # error        error message (empty string if no error occurred).
  #
  #                                first version:         Dec.  08, 2007 PR
  #                                add parameter "thin":  March 26, 2008 PR
  #                                minor modification:    March 28, 2009 PR
  # -----------------------------------------------------------------------
  
  # set up and initialize arrays:
  
  returnsize         <- floor(sampsize/thin)+1
  z.sample           <- matrix(data=NA,nrow=returnsize,ncol=length(z.ini))
  colnames(z.sample) <- names(z.ini)
  log.pdf.sample     <- rep(NA,returnsize)
  reject.freq        <- 1
  error              <- ""
  
  # calculate Cholesky decomposition of variance-covariance matrix:
  
  R <- diag(rep(1,length(prop.sd)))
  if ( is.matrix(prop.cor) ) R <- prop.cor
  if ( (nrow(R)!=length(prop.sd)) | (ncol(R)!=length(prop.sd)) )
  {
    error <- paste("markovchain.metropolis:",
                   "illegal dimension of correlation matrix")
  }
  if ( nchar(error) == 0 )
  {
    if ( length(prop.sd) == 1 ) sigma <- R * prop.sd^2
    else                        sigma <- diag(prop.sd) %*% R %*% diag(prop.sd)
    A.chol <- try(t(chol(sigma)),silent=FALSE)
    if( inherits(A.chol,"try-error") )
    {
      error <- paste("markovchain.metropolis:",
                     "unable to calculate Cholesky decomposition of variance-covariance matrix")
    }
  }
  
  # initialize Markov chain:
  
  if ( nchar(error) == 0 )
  {
    z.current         <- z.ini
    log.pdf.current   <- log.pdf(z.current,...)
    z.sample[1,]      <- z.current
    log.pdf.sample[1] <- log.pdf.current
    if ( is.na(log.pdf.sample[1]) )
    {
      error <- paste("markovchain.metropolis:",
                     "probability density is zero at initial point of chain")
    }
  }
  if ( nchar(error) == 0 )
  {
    num.accept <- 0
    for ( i in 2:returnsize )
    {
      for ( j in 1:thin )
      {
        # calculate suggested new sample point:
        
        prop.unif <- runif(length(z.ini),min=0,max=1)
        jump <- A.chol %*% qnorm(prop.unif)
        z.suggested <- z.current + as.vector(jump)
        
        # calculate log pdf at suggested new sample point:
        
        log.pdf.suggested <- log.pdf(z.suggested,...)
        
        # accept new point with probability r=pdf.suggested/pdf.prev
        
        accept <- FALSE
        if ( is.finite(log.pdf.suggested) )
        {
          if ( log.pdf.suggested > log.pdf.current )
          {
            accept <- TRUE
          }
          else
          {
            r <- exp(log.pdf.suggested-log.pdf.sample[i-1])
            if ( runif(n=1,min=0,max=1) <= r )
            {
              accept <- TRUE
            }
          }
        }
        if ( accept == TRUE )
        {
          z.current       <- z.suggested
          log.pdf.current <- log.pdf.suggested
          num.accept      <- num.accept+1
        }
        reject.freq <- ((i-2)*thin+j-num.accept)/((i-2)*thin+j)
      }
      z.sample[i,]      <- z.current
      log.pdf.sample[i] <- log.pdf.current
    }
  }
  
  # collect and return results:
  
  res <- list(z           = z.sample,
              log.pdf     = log.pdf.sample,
              reject.freq = reject.freq,
              error       = error)
  return(res)
}


################################################################################