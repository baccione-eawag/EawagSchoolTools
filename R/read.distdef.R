################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


read.distdef <- function(file,sep="\t")
{
  distdef.tab <- read.delim(file,header=F,sep=sep)
  distdef.tab <- cbind(distdef.tab,rep(NA,nrow(distdef.tab)))
  distdef <- list()
  for ( i in 1:nrow(distdef.tab) )
  {
    max.ind <- min(match(NA,t(distdef.tab[i,])),match("",t(distdef.tab[i,])),na.rm=T)-1
    distdef[[i]] <- as.vector(t(distdef.tab[i,][2:max.ind]))
  }
  names(distdef) <- as.character(distdef.tab[,1])
  return(distdef)
}

################################################################################