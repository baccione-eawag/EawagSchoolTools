################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
# last modification: Lukas M. Weber, 2014/05/14                                #
#                                                                              #
################################################################################


# sens.loc
# ========

# purpose:
# calculate the local sensitivity matrix of a model
# (matrix of partial derivatives of model results with respect to parameters)

# arguments:
# par:         named parameter vector; passed to the model as separate arguments
# model        function representing the model; typically returns a vector,
#              but could also return a matrix
# par.inc:     increments of parameters used to approximate the derivatives
# ...          further arguments are passed to model

# output:
# matrix of partial derivatives of model results (rows) 
# with respect to parameters (columns)
# or 3 dimensional array of partial derivatives if model output is a matrix

sens.loc <- function(par,model,par.inc=0.01*par,...)
{
  if ( length(par) != length(par.inc) ) 
    stop("*** error in sens.loc: unequal length of par and par.inc")
  if ( min(abs(par.inc)) == 0 )
    stop("*** error in sens.loc: elements of par.inc must not be zero")
  res.par <- model(par,...)
  V = NA
  if ( is.vector(res.par) )
  {
    V <- matrix(NA,nrow=length(res.par),ncol=length(par))
    colnames(V) <- names(par)
    rownames(V) <- names(res.par)
    for( j in 1:length(par) )
    {
      par.j <- par
      par.j[j] <- par.j[j] + par.inc[j]
      V[,j] <- ( model(par.j,...) - res.par ) / par.inc[j]
    }
  }
  else
  {
    if ( is.matrix(res.par) )
    {
      V <- array(NA,dim=c(dim(res.par),length(par)),
                 dimnames=list(rownames(res.par),
                               colnames(res.par),
                               names(res.par)))
      for ( j in 1:length(par) )
      {
        par.j <- par
        par.j[j] <- par.j[j] + par.inc[j]
        V[,,j] <- ( model(par.j,...) - res.par ) / par.inc[j]
      }
    }
  }
  return(V)
}


################################################################################


# sens.loc.explicitpars
# =====================

# Modified version of "sens.loc" that calls a model with explicitly listed 
# arguments (this is primarily for internal use when internally calling nls)
#
# model.name is the the name of the function with explicit parameter arguments

sens.loc.explicitpars <- function(par,model.name,par.inc=0.01*par,x)
{
  if ( length(par) != length(par.inc) ) 
    stop("*** error in sens.loc.explicitpars: unequal length of par and par.inc")
  if ( min(abs(par.inc)) == 0 )
    stop("*** error in sens.loc.explicitpars: elements of par.inc must not be zero")
  args <- as.list(par)
  if ( length(x) > 0 ) args <- c(args,x)
  res.par <- do.call(model.name,args=args)
  V <- matrix(NA,nrow=length(res.par),ncol=length(par))
  colnames(V) <- names(par)
  rownames(V) <- names(res.par)
  for( j in 1:length(par) )
  {
    par.j <- par
    par.j[j] <- par.j[j] + par.inc[j]
    args <- as.list(par.j)
    if ( length(x) > 0 ) args <- c(args,x)
    res.j <- do.call(model.name,args=args)
    V[,j] <- ( res.j - res.par ) / par.inc[j]
  }
  return(V)
}


################################################################################


# sens.var.samp
# =============

# function to calculate the first order variance based sensitivity
# coefficients from a parameter sample and the corresponding result sample

# Arguments:
# ----------
#
# parsamp	   matrix containing the parameter sample (each row corresponds to
#            a sample point)						
#													
# ressamp	   vector (of only one result per parameter sample point) or matrix 
#            of results corresponding to the parameter sample (each row 
#            provides the results corresponding to the parameter values in 
#            the same row of parsamp)									
#													
# nbin 	     number of quantile intervals for which the conditional means	
#		         are calculated, default is the square root of the sample size								
#
# method     "smooth", "loess", "lpepa", "lpridge"	
#            routine to be used for smoothing
#
# order		   order of local regression polynomial or of kernel
#
# bandwidth  method="lpepa" or method="lpdidge" only: bandwidth of the
#            smoothing algorithm
#
# span       method="loess" only: fration of points used for local regression 
#
# sd.rel     method="smooth" only: standard deviation of normal distribution of 
#            smoothing algorithm relative to the 99% quantile interval
#
#	plot       logical variable indicating if a scatter plot of the relationship
#            between parameter and model output should be plotted								
#													
# Output:												
# -------												
#													
# List with two elements:	
#								
# var	       vector with the total variance of each column (each model	 
#            output) of the ressamp matrix						
#													
# var.cond.E matrix with the variance of the conditional expected value	
#			       of each parameter (columns) and each model output (rows)	

sens.var.samp <- function(parsamp,ressamp,nbin=NA,
                          method="smooth",order=2,
                          bandwidth=1,span=0.75,sd.rel=0.1,
                          plot=F)
{
  #package("lpridge")
    ## packages loaded via NAMESPACE: changed LW
  
  if ( is.vector(ressamp) ) ressamp <- as.matrix(ressamp,ncol=1)
  npar  <- ncol(parsamp)
  nsamp <- nrow(parsamp)
  nres  <- ncol(ressamp)
  if ( is.na(nbin) ) nbin <- ceiling(sqrt(nsamp))
  if ( nrow(parsamp ) != nrow(ressamp)) 
  {
    stop ("ressamp and parsamp do not have the same row length")
  }
  
  var_k <- rep(NA,nres)
  names(var_k) <- colnames(ressamp)
  for ( k in 1:nres ) var_k[k] <- var(ressamp[,k])
  
  var_q <- data.frame(matrix(NA,nrow=nres,ncol=npar))
  colnames(var_q) <- colnames(parsamp)
  rownames(var_q) <- colnames(ressamp)
  
  for ( i in 1:npar )
  {
    q  <- quantile(parsamp[,i],probs=(((1:nbin)-0.5)/nbin),
                   na.rm=FALSE,names=TRUE,type=7)
    for ( k in 1:nres )
    {
      if ( method == "lpepa" )
      {
        mean_res <- lpepa(parsamp[,i],ressamp[,k],bandwidth=bandwidth,
                          x.out=q,order=order)$est
      }
      else
      {
        if ( method == "lpridge" )
        {
          mean_res <- lpridge(parsamp[,i],ressamp[,k],bandwidth=bandwidth,
                              x.out=q,order=order)$est
        }
        else
        {
          if ( method == "loess" )
          {
            #mean_res <- lowess(parsamp[,i],ressamp[,k],f,
            #                   x.out=q,korder=order,hetero=T)$est
            res <- 
              loess(y~x,
                    data=data.frame(x=parsamp[,i],y=ressamp[,k]),
                    span=span,degree=order)
            mean_res <- predict(res,newdata=data.frame(x=q))
          }
          else
          {
            if ( method == "smooth" )
            {
              if ( order == 2 ) m <- "quadratic"
              else              m <- "linear"
              mean_res <- Smooth(parsamp[,i],ressamp[,k],
                                 sigma=((q[nbin]-q[1])*sd.rel),
                                 newx=q,
                                 method=m)$y
            }
            else
            {
              stop(paste("calc.var.sens: method",
                         method,
                         "not implemented"))
            }
          }
        }
      }
      
      var_q[k,i] <- var(mean_res)
      if ( plot )
      {
        plot(parsamp[,i],ressamp[,k],pch=19,cex=0.2,
             xlab=colnames(parsamp)[i],ylab="Y")
        lines(q,mean_res,lwd=3,col="red")
      }
    }
  }
  
  return (list(var=var_k,var.cond.E=var_q))
}


################################################################################