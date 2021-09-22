################################################################################
#                                                                              #
# Functions for R package "bayesBias" : Bayesian Inference with Consideration  #
# of Bias                                                                      #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# Function implementing a simple standard posterior probability density (up to
# a multiplicative normalizing constant). It is based on a likelihood 
# definition provided by the user and a multivariate normal or
# lognormal distribution at the prior.
# This log posterior implementation serves as a template for implementing more
# specific posteriors.

Logposterior <- function(par,model,L,y,
                         prior.dist="lognormal",prior.mean=1,prior.sd=1,
                         prior.cor=NA,prior.def=NA,
                         loglikeli=Loglikeli,
                         ...)
{
  logprior  <- calcpdf_mv(z=par,dist=prior.dist,mean=prior.mean,
                          sd=prior.sd,cor=prior.cor,distdef=prior.def)
  loglik <- rep(NA,length(logprior))
  if ( length(logprior) == 1 )
  {
    if ( !is.na(logprior) )
    {
      loglik <- loglikeli(par=par,model=model,L=L,y=y,...)
    }
  }
  else
  {
    ind <- !is.na(logprior)
    if ( sum(ind) > 0 )
    {
      loglik[ind] <- loglikeli(par=par[ind,],model=model,L=L,y=y,...)
    }
  }
  return(logprior+loglik)
}


################################################################################