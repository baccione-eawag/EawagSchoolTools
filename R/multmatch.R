################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


multmatch <- function(x,table,incomparables=NULL)
{
  table.loc <- table
  inds <- numeric(0)
  while (TRUE)
  {
    ind <- match(x,table.loc,incomparables)
    if ( is.na(ind) )
    {
      if ( length(inds) == 0 ) return(ind)
      else                     return(inds)
    }
    inds <- c(inds,ind)
    table.loc[ind] <- paste(x,"_",sep="")
  }  
}


################################################################################