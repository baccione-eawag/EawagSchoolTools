################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# generate a random sample from a uniform distribution in a ball of given radius

# see introduction of Harman, R. and Lacko, V., On decompositional algorithms
# for uniform sampling from n-spheres and n-balls, Journal of Multivariate 
# Analysis 101, 2297-2304, 2010 for details and references 

randsamp.ball <- function(sampsize=1,dim=1,radius=1)
{
  res <- matrix(data=rnorm(sampsize*dim),nrow=sampsize)
  norm <- function(x) { return(sqrt(sum(x^2))) }
  rnorms <- apply(res,1,norm)
  res <- diag(radius/rnorms) %*% res
  res <- diag(runif(sampsize)^(1/dim)) %*% res
  return(res) 
}

################################################################################