################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# Confint
# =======

# purpose:
# calculate confidence intervals based on the Student t distribution
# or on the normal distribution (if df=NA)

# arguments:
# est:         vector of point estimates
# sd:          vector of estimated standard deviations
# df:          degrees of freedom
# level:       probability level of confidence intervals (default 0.95)

# output:
# table of lower and upper bounds of confidence intervals

Confint <- function(est,sd,df,level)
{
  alpha <- 1- level
  confint <- matrix(NA,nrow=length(est),ncol=2)
  colnames(confint) <- c(paste(100*alpha/2,"%",sep=""),paste(100*(1-alpha/2),"%",sep=""))
  rownames(confint) <- names(est)
  if ( is.na(df) )
  {
    fact <- qnorm(1-alpha/2)
  }
  else
  {
    fact <- qt(1-alpha/2,df)
  }
  confint[,1] <- est - fact*sd
  confint[,2] <- est + fact*sd
  return(confint)
}


################################################################################