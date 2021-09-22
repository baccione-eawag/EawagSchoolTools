################################################################################
#                                                                              #
# Functions for R package "bayesBias" : Bayesian Inference with Consideration  #
# of Bias                                                                      #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# function to improve jump distribution:

Proposal.dist <- function(postsamp,fact.sd,fract.burnin=0,fact.cor=1,plot=F)
{
  ind.end   <- nrow(postsamp)
  ind.start <- as.integer(fract.burnin*(ind.end-1))+1
  sd.postsamp <- apply(postsamp,2,sd)
  postsamp.local <- postsamp[ind.start:ind.end,sd.postsamp!=0]
  if ( is.vector(postsamp.local) )
  {
    postsamp.local <- as.matrix(postsamp.local,nrow=length(postsamp.local))
    colnames(postsamp.local) <- colnames(postsamp)[sd.postsamp!=0]
    postsamp.local <- data.frame(postsamp.local)
  }
  sd  <- fact.sd*apply(postsamp.local,2,sd)
  cor <- NA
  if ( ncol(postsamp.local) > 1 )
  {
    corr <- cor(postsamp.local)
    if ( plot )
    {
      image(x=1:ncol(corr),y=1:nrow(corr),z=abs(corr),zlim=c(0,1),
            col=grey((100:0)/100),xlab="variable index",ylab="variable index",
            main="structure of correlation matrix")
      abline(v=0.5)
      abline(h=0.5)
      abline(v=ncol(corr)+0.5)
      abline(h=nrow(corr)+0.5)
    }
    corr <- fact.cor*corr
    diag(corr) <- rep(1,nrow(corr))
    try(chol(corr)) # test if Cholesky factorization works 
    # (important for subsequent sampling)
  }
  return(list(sd=sd,cor=corr))
}


################################################################################