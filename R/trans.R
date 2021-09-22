################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# trans.to.interval
# =================

# purpose:
# transforms the real axis to an interval (default: unit interval)

# arguments:
# x:       data to be transformed
# min:     minimum of the interval (default: 0)
# max:     maximum of the interval (default: 1)

# output:
# transformed data

trans.to.interval <- function(x,min=0,max=1)
{
  y <- 0.5*(min+max) + (max-min)/pi*atan(x)
  return(y)
}


################################################################################


# trans.from.interval
# ===================

# purpose:
# transforms an interval (default: unit interval) to the real axis

# arguments:
# y:       data to be transformed
# min:     minimum of the interval (default: 0)
# max:     maximum of the interval (default: 1)

# output:
# transformed data

trans.from.interval <- function(y,min=0,max=1)
{
  x <- tan(0.5*pi*(2*y-max-min)/(max-min))
  return(x)
}


################################################################################


# trans.par.normal.tovec
# ======================

trans.par.normal.tovec <- function(mean,sd,cor,trans=T,max.cor=0.5)
{
  n <- length(mean)
  
  par <- rep(NA,n*(n+3)/2)
  par[1:n] <- mean
  par[(n+1):(2*n)] <- sd
  if ( n > 1 )
  {
    k <- 2*n
    for ( i in 1:(n-1) ) 
    {
      par[k+1:(n-i)] <- cor[i,(i+1):n]
      k <- k + n - i
    }
  }
  names(par) <- c(names(mean),names(mean),rep("cor",n*(n-1)/2))
  
  if ( trans )
  {
    par[(n+1):(2*n)] <- log(par[(n+1):(2*n)])
    par[(2*n+1):(n*(n+3)/2)] <- 
      trans.from.interval(par[(2*n+1):(n*(n+3)/2)],
                          min=-max.cor,max=max.cor)
  }
  
  return(par)
}


################################################################################


# trans.par.normal.fromvec
# ========================

trans.par.normal.fromvec <- function(par,trans=T,max.cor=0.5)
{
  n <- (-3+sqrt(9+8*length(par)))/2
  
  if ( length(par) != n*(n+3)/2 )
  {
    cat("trans.par.normal.fromvec:",
        "illegal length of parameter vector:",length(par),"\n")
    mean <- NA
    sd   <- NA
    cor  <- NA
  }
  else
  {
    if ( trans )
    {
      par[(n+1):(2*n)] <- exp(par[(n+1):(2*n)])
      par[(2*n+1):(n*(n+3)/2)] <- 
        trans.to.interval(par[(2*n+1):(n*(n+3)/2)],
                          min=-max.cor,max=max.cor)
    }
    
    mean <- par[1:n]
    sd   <- par[(n+1):(2*n)]
    cor  <- diag(rep(1,n),nrow=n)
    if ( n > 1 )
    {
      k <- 2*n
      for ( i in 1:(n-1) )
      {
        cor[i,(i+1):n] <- par[k+1:(n-i)]
        cor[(i+1):n,i] <- par[k+1:(n-i)]
        k <- k + n - i
      }
    }
    names(mean)   <- names(par)[1:n]
    names(sd)     <- names(par)[1:n]
    rownames(cor) <- names(par)[1:n]
    colnames(cor) <- names(par)[1:n]
  }
  
  return(list(mean=mean,sd=sd,cor=cor))
}            


################################################################################