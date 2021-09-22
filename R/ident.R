################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# calculation of index combinations (auxiliary function used in ident)
# ====================================================================

comb <- function(n,p)
{
  # -----------------------------------------------------------------------
  # This function calculates all combination of subsets of length p out
  # of n indices.
  #
  # Arguments:
  # n:   number of indices.
  # p:   length of subset of indices.
  #
  # Return Value:
  # matrix with subsets of length p as rows.
  #
  #                                         Peter Reichert    Dec. 27, 2002
  # -----------------------------------------------------------------------
  
  # check input:
  if ( p > n ) stop("comb: illeal arguments (p>n)")
  
  # initialize array and auxiliary variables:
  num.comb <- choose(n,p)
  comb <- matrix(nrow=num.comb,ncol=p)
  ind <- 1:p
  pointer <- p
  
  # calculate index combinations:
  for ( i in 1:num.comb )
  {
    comb[i,] <- ind
    ind[pointer] <- ind[pointer]+1
    if ( ind[pointer] > n )
    {
      while ( pointer > 1 )
      {
        pointer <- pointer-1
        if ( ind[pointer] < n-(p-pointer) )
        {
          ind[pointer] <- ind[pointer]+1
          for ( j in (pointer+1):p )
          {
            ind[j] <- ind[pointer]+j-pointer
          }
          pointer <- p
          break
        }
      }
    }
  }
  
  # return results:
  return(comb)
}


################################################################################


# calculation of collinearity index (auxiliary function used in ident)
# ====================================================================

collind <- function(sen.scaled)
{
  # -----------------------------------------------------------------------
  # This function calculates the collinearity index from a scaled 
  # sensitivity matrix.
  #
  # Arguments:
  # sen:       matrix of model sensitivities (scaled partial derivatives
  #            of model outcomes with respect to model parameters:
  #            delta.par/scale dy/dpar); the columns of sen refer to
  #            different model parameters, the rows to different model
  #            outcomes.
  #
  # Return Value:
  # collinearity index (real number).
  #
  #                                         Peter Reichert    Dec. 27, 2002
  # -----------------------------------------------------------------------
  
  # normalize sensitivity functions:
  num.par <- ncol(sen.scaled)
  norms <- numeric(num.par)
  for ( i in 1:num.par )
  {
    norms[i] <- sqrt( sen.scaled[,i] %*% sen.scaled[,i] )
  }
  sen.norm <- sen.scaled %*% diag( 1/norms )
  
  # calculate collinearity index:
  collind <- 1/sqrt(min(eigen( t(sen.norm) %*% sen.norm )$values))
  
  # return result:
  return(collind)
}


################################################################################


# calculation of identifiability measures
# =======================================

ident <- function(sen,delta.par=0,scale=0,max.subset.size=0)
{
  # -----------------------------------------------------------------------
  # This function calculates a parameter sensitivity ranking and 
  # collinearity indices for a series of parameter combinations
  # based on liniear sensitivity functions of a model, parameter
  # uncertainty ranges and scale factors of model results.
  # This function is a simplified version of the program "ident"
  # available at http://www.ident.eawag.ch
  #
  # Arguments:
  # sen:       matrix of model sensitivities (partial derivatives
  #            of model outcomes with respect to model parameters:
  #            dy/dpar); the columns of sen refer to different 
  #            model parameters, the rows to different model outcomes.
  # delta.par: model parameter uncertainty ranges (length equal to 
  #            the number of columns of sen); if zero, then all ranges
  #            are assumed to be unity.
  # scale:     scaling factors of model results (if zero, then all 
  #            scaling factors are assumed to be unity).
  #
  # Return Value:
  # List of delta.msqr, collind.
  #
  #                                         Peter Reichert    Dec. 27, 2002
  # -----------------------------------------------------------------------
  
  # determine number of model parameters:
  num.out <- nrow(sen)
  num.par <- ncol(sen)
  names.par <- colnames(sen)
  if ( length(names.par) != num.par ) names.par <- paste("par",1:num.par,sep="")
  if ( max.subset.size == 0 ) max.subset.size <- min(num.par,4)
  
  # apply parameter uncertainty ranges and scale factors if available:
  sen.scaled <- sen
  if ( length(delta.par) == num.par ) sen.scaled <- sen.scaled %*% diag(delta.par)
  if ( length(scale)     == num.out ) sen.scaled <- diag(1/scale) %*% sen.scaled
  
  # calculate sensitivity ranking:
  delta.msqr <- numeric(num.par)
  names(delta.msqr) <- names.par
  for ( i in 1:num.par )
  {
    delta.msqr[i] <- sqrt( t(sen.scaled[,i]) %*% sen.scaled[,i] ) / sqrt(num.out)
  }
  res <- list(delta.msqr=delta.msqr)
  
  if ( max.subset.size > 1 )
  {
    for ( i in 2:min(max.subset.size,num.par) )
    {
      ind <- comb(num.par,i)
      collind <- numeric(nrow(ind))
      par.set <- matrix(nrow=nrow(ind),ncol=i)
      colnames(par.set) <- paste("par.",1:i,sep="")
      for ( j in 1:nrow(ind) )
      {
        collind[j] <- collind(sen.scaled[,ind[j,]])
        for ( k in 1:i )
        {
          par.set[j,k] <- names.par[ind[j,k]]
        }
      }
      if ( nrow(par.set) > 1 )
      {
        ind.sorted <- order(collind)
        res[[paste("collind.",i,sep="")]] <- 
          data.frame(par.set[ind.sorted,],collind=collind[ind.sorted])
      }
      else
      {
        res[[paste("collind.",i,sep="")]] <- 
          data.frame(par.set,collind=collind)
      }
    }
  }
  
  # return results:
  return(res)
}


################################################################################