################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# Hessian
# =======

Hessian <- function(fn,par,par.inc=NA,...)
{
  if ( is.na(par.inc[1]) ) par.inc <- 0.01*par
  n <- length(par)
  h <- matrix(NA,nrow=n,ncol=n)
  colnames(h) <- names(par)
  rownames(h) <- names(par)
  trace <- matrix(NA,nrow=2*n^2+1,ncol=n+1)
  colnames(trace) <- c(names(par),"fn")
  minimum <- TRUE
  maximum <- TRUE
  par.0 <- par
  res.0 <- fn(par.0,...)
  num.eval <- 1
  trace[1,] <- c(par.0,res.0)
  row.trace <- 1
  for ( i in 1:n )
  {
    for ( j in 1:i )
    {
      if ( i==j )
      {
        par.inc.tmp <- par.inc
        counter <- 0
        while (TRUE)
        {
          par.1 <- par.0
          par.1[i] <- par.1[i] + par.inc.tmp[i]
          res.1 <- fn(par.1,...)
          num.eval <- num.eval+1
          if ( !is.na(res.1) )
          { 
            if ( res.1 > res.0 ) maximum <- FALSE
            if ( res.1 < res.0 ) minimum <- FALSE
            par.2 <- par.0
            par.2[i] <- par.2[i] - par.inc.tmp[i]
            res.2 <- fn(par.2,...)
            num.eval <- num.eval+1
            if ( !is.na(res.2) )
            {
              if ( res.2 > res.0 ) maximum <- FALSE
              if ( res.2 < res.0 ) minimum <- FALSE
              h[i,i] <- (res.1 - 2*res.0 + res.2)/par.inc.tmp[i]^2
              trace[row.trace+1,] <- c(par.1,res.1)
              trace[row.trace+2,] <- c(par.2,res.2)
              row.trace <- row.trace + 2
              break
            }
          }
          counter <- counter + 1
          if ( counter > 5 ) stop("hessian: unable to calculate hessian")
          par.inc.tmp <- 0.5*par.inc.tmp
        }
      }
      else
      {
        par.inc.tmp <- par.inc
        counter <- 0
        while (TRUE)
        {
          par.1 <- par.0
          par.1[i] <- par.1[i] + par.inc.tmp[i]
          par.1[j] <- par.1[j] + par.inc.tmp[j]
          res.1 <- fn(par.1,...)
          num.eval <- num.eval + 1
          if ( !is.na(res.1) )
          {
            if ( res.1 > res.0 ) maximum <- FALSE
            if ( res.1 < res.0 ) minimum <- FALSE
            par.2 <- par.0
            par.2[i] <- par.2[i] + par.inc.tmp[i]
            par.2[j] <- par.2[j] - par.inc.tmp[j]
            res.2 <- fn(par.2,...)
            num.eval <- num.eval + 1
            if ( !is.na(res.2) )
            {
              if ( res.2 > res.0 ) maximum <- FALSE
              if ( res.2 < res.0 ) minimum <- FALSE
              par.3 <- par.0
              par.3[i] <- par.3[i] - par.inc.tmp[i]
              par.3[j] <- par.3[j] + par.inc.tmp[j]
              res.3 <- fn(par.3,...)
              num.eval <- num.eval + 1
              if ( !is.na(res.3) )
              {
                if ( res.3 > res.0 ) maximum <- FALSE
                if ( res.3 < res.0 ) minimum <- FALSE
                par.4 <- par.0
                par.4[i] <- par.4[i] - par.inc.tmp[i]
                par.4[j] <- par.4[j] - par.inc.tmp[j]
                res.4 <- fn(par.4,...)
                num.eval <- num.eval + 1
                if ( !is.na(res.4) )
                {
                  if ( res.4 > res.0 ) maximum <- FALSE
                  if ( res.4 < res.0 ) minimum <- FALSE
                  h[i,j] <- (res.1 - res.2 - res.3 + res.4)/
                    (4*par.inc.tmp[i]*par.inc.tmp[j])
                  h[j,i] <- h[i,j]
                  trace[row.trace+1,] <- c(par.1,res.1)
                  trace[row.trace+2,] <- c(par.2,res.2)
                  trace[row.trace+3,] <- c(par.3,res.3)
                  trace[row.trace+4,] <- c(par.4,res.4)
                  row.trace <- row.trace + 4
                  break
                }
              }
            }
          }
          counter <- counter + 1
          if ( counter > 5 ) stop("hessian: unable to calculate hessian")
          par.inc.tmp <- 0.5*par.inc.tmp
        }
      } 
    }
  }
  return(list(h        = h,
              minimum  = minimum,
              maximum  = maximum,
              num.eval = num.eval,
              trace    = trace))
}      


################################################################################