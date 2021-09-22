################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
# last modification: Lukas M. Weber, 2014/05/14                                #
#                                                                              #
################################################################################


# Plot.res
# ========

# purpose:
# plot results provided as a numeric vector with result codes

# arguments:
# res:         vector of results
# L:           observation layout
# xlim:        optional limits of the x axis
# ylim:        optional limits of the y axis
# markers:     if TRUE plot markers instead of lines
# header:      optional header of the plot
# xlab:        optional label of the x axis
# ylab:        optional label of the y axis
# pos:         position of legend (only if more than one variable encoded)

# output:
# plot of all variables as a function of the independent variable

Plot.res <- function(res,L,xlim=NA,ylim=NA,
                     markers=F,header="",
                     xlab="",ylab="",pos="topright")  
                       ## changed LW: added argument L
{
  varnames <- unique(L[,1])
  if ( is.na(xlim[1]) )      xlim <- range(L[,2])
  if ( is.na(ylim[1]) )      ylim <- range(res)
  if ( nchar(ylab[1]) == 0 ) ylab <- paste(varnames,collapse=", ")
  plot(numeric(0),numeric(0),type="n",xlim=xlim,ylim=ylim,
       xlab=xlab,ylab=ylab,main=header)
  if ( markers )
  {
    for ( i in 1:length(varnames) )
    {
      ind <- L[,1] == varnames[i]
      points(L[ind,1],res[ind],pch=i)
    }
    if ( length(varnames) > 1 )
    {
      legend(x=pos,legend=varnames,pch=1:length(varnames))
    }
  }
  else
  {
    for ( i in 1:length(varnames) )
    {
      ind <- L[,1] == varnames[i]
      lines(L[ind,2],res[ind],lty=i)
    }
    if ( length(varnames) > 1 )
    {
      legend(x=pos,legend=varnames,lty=1:length(varnames))
    }
  }
}


################################################################################


# Plot.vars
# =========

# purpose:
# plot variables provided as a data frame or matrix, with the corresponding 
# variable information provided by an observation layout

# arguments:
# vars:        matrix or data frame with variables, and variable information 
#              provided by the observation layout
# L:           observation layout
# ncol:        optional number of columns of sub-panels of the plot
# mar:         optional specification of margins in the form 
#              c(bottom,left,top,right)
# ylim:        optional named (by variable name) list of limits of the y axes
# markers:     if TRUE plot markers instead of lines
# header:      optional named (by variable name) list of headers of the plots
# xlab:        optional label of the x axis
# ylab:        optional label of the y axis
# pos:         position of legend (only if more than one variable)

# output:
# plot of all variables as a function of the independent variable
# (note that variable names and values of the independent variable are
# encoded in the names of the components of the result vector)

Plot.vars <- function(vars,L,ncol=NA,mar=NA,
                      ylim=list(),markers=F,
                      headers=list(),xlab="",ylab="",pos="topright")
                        ## changed LW: added argument L
{
  nvar <- ncol(vars)
  if ( is.na(ncol) ) nc <- floor(sqrt(nvar))
  nr <- ceiling(nvar/nc)
  marg <- mar
  if ( is.na(marg[1]) ) marg <- c(4.5,4.0,2.5,1.0) # c(bottom, left, top, right)
  par.def <- par(no.readonly=T)
  par(mfrow=c(nr,nc),xaxs="i",yaxs="i",mar=marg) 
  for ( i in 1:nvar )
  {
    name <- colnames(vars)[i]
    ylim.i <- c(min(vars[,i]),max(vars[,i]))
    if ( ylim.i[1] == ylim.i[2] ) 
    {
      ylim.i[1] <- 0.5*ylim.i[1]
      ylim.i[2] <- 1.5*ylim.i[2]
    }
    if ( length(ylim) > 0 )
    {
      ind <- match(name,names(ylim))
      if ( !is.na(ind) )
      {
        ylim.i <- ylim[[ind]]
      }
    }
    header.i <- name
    if ( length(headers) > 0 )
    {
      ind <- match(name,names(headers))
      {
        if ( !is.na(ind) )
        {
          header.i <- headers[[ind]]
        }
      }
    }
    res <- as.numeric(vars[,i])
    names(res) <- rownames(vars)
    Plot.res(res,
             L,
             header=header.i,markers=markers,
             xlab=xlab,ylab=ylab,pos=pos,
             xlim=NA,ylim=ylim.i)
  }
  par(par.def)
}


################################################################################


# Contour of the probability density function of a bivariate normal or 
# lognormal distribution

contourpdf <- function(calcpdf_mv,norm=T,xlim=c(-3,3),ylim=c(-3,3),
                       levels=c(0.05,0.5,0.95),res=20,lty="solid",
                       ...)
{
  # -----------------------------------------------------------------------
  # This function plots contour lines of a normalized or unnormalized
  # bivariate probability density function. If the function is not 
  # normalized, the integral over the domain specified by xlim and ylim is  
  # used for normalization (this leads to incorrect results if this domain 
  # does not contain most of the distribution).
  #
  # Arguments:
  # calcpdf_mv: function to calculate the log pdf values at a set of
  #             locations specified by its first argument. Further 
  #             arguments specified under ... will be passed to this 
  #             function
  # norm:       TRUE if the probability density is normalized,
  #             FALSE if integration over the domain given by xlim and ylim
  #             should be used for normalization
  # xlim:       bounds of the integration range for the first variable
  # ylim:       bounds of the integration range for the second variable 
  # levels:     vector of probabilities to be contained in the contour line
  # res:        resolution of grid used to countour
  #             (number of points in each dimension)
  # lty:        line type of contour lines
  #
  # Return Value:
  # integral of the probability density over the given range at the given
  # resolution
  #
  #                                        Peter Reichert    Feb.  02, 2005
  #                                        last modification March 31, 2008
  # -----------------------------------------------------------------------
  
  dx.grid <- (xlim[2]-xlim[1])/res
  dy.grid <- (ylim[2]-ylim[1])/res
  x.grid <- seq(xlim[1]+dx.grid/2,xlim[2]-dx.grid/2,len=res)
  y.grid <- seq(ylim[1]+dy.grid/2,ylim[2]-dy.grid/2,len=res)
  xarray <- vector()
  for ( i in 1:res ) xarray = c(xarray,rep(x.grid[i],res))
  yarray <- rep(y.grid,res)
  z.sample <- cbind(xarray,yarray)
  logpdf <- calcpdf_mv(z.sample,...)
  if ( norm == F) logpdf <- logpdf-max(logpdf,na.rm=T)
  pdf <- ifelse(is.na(logpdf),0,exp(logpdf))
  pdf.sort <- sort(pdf,decreasing=T)
  integral.cum <- pdf.sort*dx.grid*dy.grid
  for ( i in 2:(res*res) ) 
  {
    integral.cum[i] <- integral.cum[i-1]+integral.cum[i]
  }
  integral <- integral.cum[res*res]
  if ( norm == F) integral.cum <- integral.cum/integral
  logpdflevels <- vector()
  index <- 1
  levelssort <- sort(levels)
  for ( i in 1:(res*res) )
  {
    if ( integral.cum[i] > levelssort[index] )
    {
      logpdflevels[index] <- log(pdf.sort[i])
      index <- index + 1
      if ( index > length(levelssort) ) break
    }
  }
  logpdf <- matrix(logpdf,nrow=res,ncol=res,byrow=TRUE)
  contour(x=x.grid,y=y.grid,z=logpdf,levels=logpdflevels,
          add=TRUE,drawlabels=FALSE,lty=lty)
  return(integral)
}


################################################################################


# function to plot Markov chains:
# -------------------------------

Plot.chains <- function(postsamp,ncol=NA,mar=NA,
                        ylim=list(),
                        titles=list(),xlab="chain index",ylab=list())
{
  nvar <- ncol(postsamp)
  if ( is.na(ncol) ) nc <- floor(sqrt(nvar))
  nr <- ceiling(nvar/nc)
  marg <- mar
  if ( is.na(marg[1]) ) marg <- c(2.0,2.0,2.5,0.5) # c(bottom, left, top, right)
  par.def <- par(no.readonly=T)
  par(mfrow=c(nr,nc),xaxs="i",yaxs="i",mar=marg) 
  for ( i in 1:nvar )
  {
    name <- colnames(postsamp)[i]
    data <- postsamp[,i]
    ylim.i <- c(min(data),max(data))
    if ( ylim.i[1] == ylim.i[2] ) 
    {
      ylim.i[1] <- 0.5*ylim.i[1]
      ylim.i[2] <- 1.5*ylim.i[2]
    }
    if ( length(ylim) > 0 )
    {
      ind <- match(name,names(ylim))
      if ( !is.na(ind) )
      {
        ylim.i <- ylim[[ind]]
      }
    }
    title.i <- name
    if ( length(titles) > 0 )
    {
      ind <- match(name,names(titles))
      {
        if ( !is.na(ind) )
        {
          title.i <- titles[[ind]]
        }
      }
    }
    ylab.i <- name
    if ( length(ylab) > 0 )
    {
      ind <- match(name,names(ylab))
      if ( !is.na(ind) )
      {
        ylab.i <- ylab[[ind]]
      }
    }
    plot(data,ylim=ylim.i,type="l",main=title.i,xlab=xlab,ylab=ylab.i)
  }
  par(par.def)
}


################################################################################


# function to plot prior and posterior marginals:
# -----------------------------------------------

Plot.margs <- function(postsamp,pridist=list(),val=NA,ncol=NA,nrow=NA,mar=NA,
                       xlim=list(),ymax=list(),
                       titles=list(),xlab=list(),ylab=list(),
                       adjust=1,lty=NA,col.post=NA,col.pri=NA)
{
  # transform samples to a list of all samples to plot the marginals of:
  
  postsamp.list <- list()
  if ( is.data.frame(postsamp) )
  {
    postsamp.list[[1]] <- postsamp
  }
  else
  {
    if ( is.matrix(postsamp) )
    {
      postsamp.list[[1]] <- postsamp
    }
    else
    {
      if ( is.list(postsamp) )  # be careful: a data frame is a list!
      {
        postsamp.list <- postsamp
      }
      else
      {
        stop("Plot.margs: postsamp is of illegal type")
      }
    }
  }
  nsamp <- length(postsamp.list)
  
  # transform prior definitions to a list of priors to plot:
  
  pridist.list <- list()
  if ( length(pridist) > 0 )
  {
    if ( length(names(pridist)) != length(pridist) )
    {
      pridist.list <- pridist
    }
    else
    {
      pridist.list[[1]] <- pridist
    }
  }
  npri <- length(pridist.list)
  
  # get all variable names of the sample: 
  
  var <- colnames(postsamp.list[[1]])
  if ( nsamp > 1 )
  {
    for ( j in 2:nsamp ) 
    {
      var <- c(var,colnames(postsamp.list[[j]]))
    }
    var <- unique(var)
  }
  nvar <- length(var)
  
  # define layout of plot panels:
  
  if ( !is.na(ncol) ) nc <- ncol
  else                nc <- floor(sqrt(nvar))
  if ( !is.na(nrow) ) nr <- nrow
  else                nr <- ceiling(nvar/nc)
  marg <- mar
  if ( is.na(marg[1]) ) marg <- c(2.0,2.0,2.5,0.5) # c(bottom, left, top, right)
  
  # define line types:
  
  lty.loc <- lty
  if ( is.na(lty.loc[1]) ) lty.loc <- 1
  if ( length(lty.loc) < nsamp+npri ) lty.loc <- 1:(nsamp+npri)
  
  # plot marginals of samples and priors:
  
  par.def <- par(no.readonly=T)
  par(mfrow=c(nr,nc),xaxs="i",yaxs="i",mar=marg) # c(bottom, left, top, right)
  for ( i in 1:nvar )
  {
    name <- var[i]
    data.min    <- NA
    data.max    <- NA
    density.max <- NA
    marg.samp <- as.list(rep(NA,nsamp))
    for ( j in 1:nsamp )
    {
      ind <- match(var[i],colnames(postsamp.list[[j]]))
      if ( !is.na(ind) )
      {
        data.min       <- min(c(data.min,postsamp.list[[j]][,ind]),na.rm=TRUE)
        data.max       <- max(c(data.max,postsamp.list[[j]][,ind]),na.rm=TRUE)
        marg.samp[[j]] <- density(postsamp.list[[j]][,ind],adjust=adjust)
        density.max    <- max(c(density.max,marg.samp[[j]]$y),na.rm=TRUE)
      }
    }
    
    xlim.i <- c(data.min,data.max)
    if ( xlim.i[1] == xlim.i[2] ) 
    {
      xlim.i[1] <- 0.5*xlim.i[1]
      xlim.i[2] <- 1.5*xlim.i[2]
    }
    if ( length(xlim) > 0 )
    {
      ind <- match(name,names(xlim))
      if ( !is.na(ind) )
      {
        xlim.i <- xlim[[ind]]
      }
    }
    
    if ( is.na(density.max) ) ylim.i <- c(0,1)
    else                      ylim.i <- c(0,1.1*density.max)
    if ( length(ymax) > 0 )
    {
      ind <- match(name,names(ymax))
      if ( !is.na(ind) )
      {
        ylim.i <- c(0,ymax[[ind]])
      }
    }
    
    title.i <- name
    if ( length(titles) > 0 )
    {
      ind <- match(name,names(titles))
      {
        if ( !is.na(ind) )
        {
          title.i <- titles[[ind]]
        }
      }
    }
    
    xlab.i <- name
    if ( length(xlab) > 0 )
    {
      ind <- match(name,names(xlab))
      if ( !is.na(ind) )
      {
        xlab.i <- xlab[[ind]]
      }
    }
    
    ylab.i <- "f"
    if ( length(ylab) > 0 )
    {
      ind <- match(name,names(ylab))
      if ( !is.na(ind) )
      {
        ylab.i <- ylab[[ind]]
      }
    }
    
    marg.pri <- list()
    if ( npri > 0 )
    {
      x <- seq(xlim.i[1],xlim.i[2],length=101)
      y <- numeric(0)
      for ( j in 1:npri )
      {
        marg.pri[[j]] <- list()
        if ( length(pridist.list[[j]]) > 0 )
        {
          ind <- multmatch(name,names(pridist.list[[j]]))
          if ( !is.na(ind[1]) )
          {
            for ( k in 1:length(ind) )
            {
              y <- calcpdf(x,pridist.list[[j]][[ind[k]]])
            }
          }
        }
        marg.pri[[j]]$x <- x
        marg.pri[[j]]$y <- y
      }
    }
    
    # plot marginals of single variable:
    
    # plot frame:
    
    plot(numeric(0),numeric(0),main=title.i,xlab=xlab.i,ylab=ylab.i,
         xlim=xlim.i,ylim=ylim.i,type="n")
    
    # plot areas:
    
    if ( npri > 0 & !is.na(col.pri) )
    {
      for ( j in 1:npri )
      {
        if ( !is.na(marg.pri[[j]])[[1]] )
        {
          n <- length(marg.pri[[j]]$x)
          polygon(c(marg.pri[[j]]$x[1],marg.pri[[j]]$x,marg.pri[[j]]$x[n],marg.pri[[j]]$x[1]),
                  c(0,marg.pri[[j]]$y,0,0),
                  border=NA,col=col.pri)
        }
      }
    }
    if ( !is.na(col.post) )
    {
      for ( j in 1:nsamp )
      {
        if ( !is.na(marg.samp[[j]])[[1]] )
        {
          n <- length(marg.samp[[j]]$x)
          polygon(c(marg.samp[[j]]$x[1],marg.samp[[j]]$x,marg.samp[[j]]$x[n],marg.samp[[j]]$x[1]),
                  c(0,marg.samp[[j]]$y,0,0),
                  border=NA,col=col.post)
        }
      }
    }
    
    # plot lines:
    
    abline(h=ylim.i[1])
    abline(h=ylim.i[2])
    abline(v=xlim.i[1])
    abline(v=xlim.i[2])
    if ( !is.na(val[var[i]]) ) abline(v=val[var[i]])
    if ( npri > 0 )
    {
      for ( j in 1:npri )
      {
        if ( !is.na(marg.pri[[j]])[[1]] )
        {
          lines(marg.pri[[j]],lty=lty.loc[nsamp+j])
        }
      }
    }
    for ( j in 1:nsamp )
    {
      if ( !is.na(marg.samp[[j]])[[1]] )
      {
        lines(marg.samp[[j]],lty=lty.loc[j])
      }
    }
  }
  par(par.def)
}


################################################################################