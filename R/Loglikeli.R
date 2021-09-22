################################################################################
#                                                                              #
# Functions for R package "bayesBias" : Bayesian Inference with Consideration  #
# of Bias                                                                      #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# Function implementing a simple standard likelihood function for a model
# provided by a deterministic function. This function provides the density of
# a multivariate normal distribution centered at the results of the 
# deterministic model.
# If there exists a component of the parameter vector labelled "sd_rel" 
# instead of the provided standard deviations, "error.sd" the standard 
# deviations "error.sd * par["sd_rel"]" are used. This makes it possible to
# estimate a common factor of a given variance-covariance structure.
# This likelihood implementation serves as a template for implementing more
# specific likelihood functions. It is called by "Logposterior".

Loglikeli <- function(par,model,L,y,
                      error.sd=NA,error.cor.inv=NA,...)
{
  sd.rel <- 1; if ( !is.na(par["sd_rel"]) ) sd.rel <- par["sd_rel"]
  res <- numeric(0)
  if ( !is.matrix(par) )  # single model evaluation at parameter vector par
  {
    mean <- model(par,L=L,...)
    sd <- error.sd; if ( is.na(sd[1]) ) sd <- rep(1,length(mean))
    if ( is.na(error.cor.inv[1]) )
    {
      res <- 0;
      for ( i in 1:length(mean) )
      {
        res <- res + dnorm(y[i],mean[i],sd[i]*sd.rel,log=TRUE)
      }
    }
    else 
    {
      res <- calcpdf_mv(z=y,dist="normal",
                        mean=mean,sd=sd*sd.rel,
                        cor.inv=error.cor.inv)
    }
  }
  else                    # multiple model evaluations for parameter sample
  {
    for ( i in 1:nrow(par) )
    {
      mean <- model(par[i,],L=L,...)
      sd <- error.sd; if ( is.na(sd[1]) ) sd <- rep(1,length(mean))
      if ( is.na(error.cor.inv[1]) )
      {
        res[i] <- 0;
        for ( i in 1:length(mean) )
        {
          res[i] <- res[i] + dnorm(y[i],mean[i],sd[i]*sd.rel,log=TRUE)
        }
      }
      else 
      {
        res[i] <- calcpdf_mv(z=y,dist="normal",
                             mean=mean,sd=sd*sd.rel,
                             cor.inv=error.cor.inv)
      }
    }
  }
  return(res)
}


################################################################################