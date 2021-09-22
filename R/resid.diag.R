################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# resid.diag
# ==========

# purpose:
# plot residual diagnostics plots

# arguments:
# obs:         vector of observations
# calc:        vector of calculated results

# output:
# diagnostic plots

resid.diag <- function(obs,calc,header="")
{
  # calculate range of measurements and normalized residuals:
  # ---------------------------------------------------------
  
  obs.min        <- min(obs)
  obs.max        <- max(obs)
  calc.min       <- min(calc)
  calc.max       <- max(calc)
  resid.norm     <- (obs-calc)/sd(obs-calc)
  resid.norm.abs <- abs(resid.norm)
  resid.max      <- max(abs(resid.norm.abs))
  resid.lim      <- 1.1*resid.max*c(-1,1)
  marker         <- 19
  
  # divide plot area into four panels:
  # ----------------------------------
  
  par.def <- par(no.readonly=TRUE)
  par(mfrow=c(2,2),xaxs="i",yaxs="i",mar=c(4.5,4,3,1),oma=c(0,0,2,0))
  
  # plot sequence of residuals:
  # ---------------------------
  
  plot(resid.norm,main="Sequence of Residuals",
       ylim=resid.lim,pch=marker,cex=0.8,
       xlab="Data Points",ylab="Normalized Residuals")
  lines(c(0,length(resid.norm)),c(0,0))
  
  # plot residuals as function of predicted values:
  # -----------------------------------------------
  
  plot(calc,resid.norm,main="Residuals vs. Predicted",
       xlim=c(calc.min,calc.max),ylim=resid.lim,
       xlab="Predicted",ylab="Normalized Residuals",pch=marker,cex=0.8)
  lines(c(calc.min,calc.max),c(0,0))
  res.lm <- lm(resid.norm.abs ~ calc)
  x.new <- c(calc.min,calc.max)
  y.new <- predict(res.lm,newdata=data.frame(calc=x.new))
  lines(x.new, y.new)
  lines(x.new,-y.new)
  
  # plot histogram of residuals:
  # ----------------------------
  
  hist(resid.norm,freq=FALSE,main="Hist. of Residuals",
       xlim=resid.lim,
       xlab="Normalized Residuals",ylab="Density")
  lines(seq(-3,3,by=0.1),dnorm(seq(-3,3,by=0.1)))
  lines(resid.lim,c(0,0))
  
  # normal quantile plot:
  # ---------------------
  
  lim <- max(resid.max,qnorm(1-0.5/length(obs))+0.1)
  qqnorm(resid.norm,main="Sample vs. Normal Quant.",
         xlab="Normal Quantiles",ylab="Sample Quantiles",
         pch=marker,cex=0.8,
         xlim=1.1*lim*c(-1,1),ylim=1.1*lim*c(-1,1))
  lines(c(-10,10),c(-10,10))
  
  # reset plot attributes:
  # ----------------------
  
  mtext(header,side=3,outer=T,adj=0.5,cex=1.2)
  
  par(par.def)
}


################################################################################


# resid.diag.boxcox
# =================

# purpose:
# plot residual diagnostics plots for given Box-Cox transformation parameters

# arguments:
# obs:         vector of observations
# calc:        vector of calculated results
# lambda1:     value of parameter lambda1 (exponent) of Box-Cox transformation
# lambda2:     value of parameter lambda2 (offset) of Box-Cox transformation

# output:
# diagnostic plots

resid.diag.boxcox <- function(obs,calc,lambda1=1,lambda2=1)
{
  resid.diag(BoxCox(obs,lambda1,lambda2),
             BoxCox(calc,lambda1,lambda2))
}


################################################################################