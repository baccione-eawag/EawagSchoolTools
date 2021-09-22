################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# didactical emulator (preliminary version, added May 30, 2012):
# ==============================================================

emulator.evaluate <- function(x, ...) UseMethod("evaluate")


################################################################################


emulator.create <- function(inp.design,
                            res.design,
                            par,
                            sd,
                            lambda,
                            alpha    = 2,
                            sd.smooth = 0,
                            pri.mean = emulator.pri.mean.lin,
                            pri.var  = emulator.pri.var.normal)
{
  #package("fpc")
  
  n.design <- length(res.design)
  dim.inp  <- length(lambda)
  if ( is.vector(inp.design) ) 
  {
    inp.design <- matrix(inp.design,nrow=n.design)
  }
  if ( nrow(inp.design) != n.design ) stop("number of outputs must be equal to number of inputs")
  if ( ncol(inp.design) != dim.inp )  stop("dimension of input must be equal to length of lambda")
  
  emu <- list()
  emu$inp.design  <- inp.design
  emu$res.design  <- res.design
  emu$par         <- par
  emu$sd          <- sd
  emu$lambda      <- lambda
  emu$alpha       <- alpha
  emu$sd.smooth   <- sd.smooth
  emu$pri.mean    <- pri.mean
  emu$pri.var     <- pri.var
  emu$pri.var.num <- pri.var(inp=inp.design,sd=sd,lambda=lambda,alpha=alpha,sd.smooth=sd.smooth)
  emu$delta       <- res.design - pri.mean(inp=inp.design,par=par)
  emu$inv         <- solve(emu$pri.var.num)
  #emu$inv         <- solvecov(emu$pri.var.num)$inv
  emu$v           <- emu$inv %*% emu$delta
  class(emu) <- "emulator.GASP"
  
  return(emu)
}


################################################################################


emulator.pri.mean.lin <- function(inp,par)
{
  if ( is.vector(inp) ) inp <- matrix(inp,ncol=length(par)-1)
  par <- as.numeric(par)
  
  n.inp   <- nrow(inp)
  dim.inp <- ncol(inp)
  dim.par <- length(par)
  
  if ( dim.par != dim.inp+1 ) stop("number of parameters must be dim.inp+1")
  
  res <- par[1] + inp %*% par[-1]
  rownames(res) <- 1:n.inp
  colnames(res) <- "y"
  
  return(res)
}


################################################################################


emulator.pri.var.normal <- function(inp,sd,lambda,alpha=2,sd.smooth=0,rows=NA,cols=NA)
{
  if ( is.vector(inp) ) inp <- matrix(inp,ncol=length(lambda))
  n.inp   <- nrow(inp)
  dim.inp <- ncol(inp)
  if ( length(lambda) != dim.inp ) stop("length of lambda must be dim.inp")
  if ( is.na(rows[1]) ) rows <- 1:n.inp
  if ( is.na(cols[1]) ) cols <- 1:n.inp
  var <- matrix(0,nrow=length(rows),ncol=length(cols))
  
  s <- matrix(0,nrow=length(rows),ncol=length(cols))
  for ( k in 1:dim.inp )
  {
    s <- s + (abs(inp[rows,k]%o%rep(1,length(cols))-rep(1,length(rows))%o%inp[cols,k])/
                lambda[k])^alpha
  }
  var <- sd*sd*exp(-s)
  var <- var+diag(sd.smooth*sd.smooth,nrow=n.inp)[rows,cols]
  
  colnames(var) <- cols
  rownames(var) <- rows
  return(var)
}


################################################################################


evaluate.emulator.GASP <- function(emulator,inp)
{
  if ( is.vector(inp) ) inp <- matrix(inp,ncol=length(emulator$lambda))
  n.inp <- nrow(inp)
  n.design <- nrow(emulator$inp.design)
  y   <- rep(NA,n.inp)
  var <- rep(NA,n.inp)
  for ( i in 1:n.inp )
  {
    pri.mean <- emulator$pri.mean(inp      = matrix(inp[i,],nrow=1),par=emulator$par)
    k        <- emulator$pri.var(inp       = rbind(emulator$inp.design,matrix(inp[i,],nrow=1)),
                                 sd        = emulator$sd,
                                 lambda    = emulator$lambda,
                                 alpha     = emulator$alpha,
                                 sd.smooth = emulator$sd.smooth,
                                 rows      = n.design+1,
                                 cols      = 1:n.design)
    K        <- emulator$pri.var(inp       = matrix(inp[i,],nrow=1),
                                 sd        = emulator$sd,
                                 lambda    = emulator$lambda,
                                 alpha     = emulator$alpha,
                                 sd.smooth = emulator$sd.smooth)
    y[i]   <- pri.mean + k %*% emulator$v
    var[i] <- K - k %*% emulator$inv %*% t(k)
  }
  return(list(inp=inp,y=y,var=var))
}


################################################################################