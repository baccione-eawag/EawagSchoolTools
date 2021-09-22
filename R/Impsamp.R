################################################################################
#                                                                              #
# Functions for R package "bayesBias" : Bayesian Inference with Consideration  #
# of Bias                                                                      #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# Calculate an importance sample of a distribution

Impsamp <- function(log.pdf,z,z.log.pdf,...)
{
  w <- rep(NA,nrow(z))
  log.pdf.values <- w
  for ( i in 1:nrow(z) )
  {
    par <- z[i,]
    names(par) <- colnames(z)
    log.pdf.values[i] <- log.pdf(par,...)
  }
  log.pdf.max <- max(log.pdf.values,na.rm=TRUE)
  w <- exp(log.pdf.values-log.pdf.max-z.log.pdf)
  w <- ifelse ( is.na(w),0,w )
  w <- w/sum(w)
  ess <- sum(w)^2/sum(w^2)
  return(list(z=z,w=w,ess=ess,
              log.pdf.dist=log.pdf.values,log.pdf.samp=z.log.pdf))
}


################################################################################