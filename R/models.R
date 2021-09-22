################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
# modified by Lukas M. Weber, 2014/05/14                                       #
# last modification: Peter REichert 2014/06/14                                 #
#                                                                              #
################################################################################


model.monod <- function(par,L)
{
  # model:
  # ------
  #
  #   growth rate as a function of substrate concentration:
  #
  #   r = r_max*C / (K + C)
  #
  # parameters:
  # -----------
  #
  #   r_max    maximum growth rate
  #   K        half-saturation concentration
  #
  # layout:
  # -------
  #
  #   L        observation layout containing values of concentration C at which 
  #            the growth rate is to be calculated
  #
  #            L should be a data frame with two columns. The first column 
  #            contains the names of the variables to be evaluated, and the 
  #            second column contains the time points (or other evaluation 
  #            steps) for the model. See help file for function "checkLayout" 
  #            in "bayesBias" package for more details. For example:
  #
  #              Monod model:           Growth model:
  #                  var  |   C            var  |   t
  #                ---------------       ---------------
  #                   r   |  0.1           C_M  |   1
  #                   r   |  0.2           C_M  |   2
  #                   r   |  0.3           C_M  |   3
  #                  ...  |  ...           ...  |  ...
  #                                        C_S  |   1
  #                                        C_S  |   2
  #                                        C_S  |   3
  #                                        ...  |  ...
  
  # check input:
  # ------------
  
  if ( length(par) < 2 )
  {
    stop("error in model.monod: insufficent number of parameters provided")
  }
  
  # process layout:
  # ---------------
  
  n <- nrow(L)
  var <- unique(L[,1])
  if ( length(var) > 1 | var != "r" )
  {
    stop(paste("*** error in model.monod:",
               "can only calculate results for variable \"r\";",
               "layout requests",
               paste(var,collapse=",")))
  }
  C <- L[,2]
  
  # assign parameter values:
  # ------------------------
  
  r_max <- par["r_max"]; if (is.na(r_max)) r_max <- par[1]
  K     <- par["K"];     if (is.na(K))     K     <- par[2]
  
  # calculate results:
  # ------------------
  
  r <- r_max*C/(K+C)
  
  # return results:
  # ---------------
  
  res <- as.numeric(r)
  names(res) <- paste("r_",C,sep="") 
  return(res)
}


################################################################################


model.monod.external <- function(par,L)
{
  # check input:
  # ------------
  
  if ( length(par) < 2 )
  {
    stop("error in model.monod: insufficent number of parameters provided")
  }
  
  # process layout:
  # ---------------
  
  n <- nrow(L)
  var <- unique(L[,1])
  if ( length(var) > 1 | var != "r" )
  {
    stop(paste("*** error in model.monod:",
               "can only calculate results for variable \"r\";",
               "layout requests",
               paste(var,collapse=",")))
  }
  C <- L[,2]
  
  # assign parameter values:
  # ------------------------
  
  r_max <- par["r_max"]; if (is.na(r_max)) r_max <- par[1]
  K     <- par["K"];     if (is.na(K))     K     <- par[2]
  
  # write parameter set and input to files:
  
  write.table(data.frame(par),file="monod.par.tmp",col.names=F,quote=F,sep="\t")
  write.table(data.frame(C),file="monod.input.tmp",col.names=F,row.names=F)
  
  # execute external simulation:
  
  system("model_monod monod.par.tmp monod.input.tmp monod.out.tmp")
  
  # read results from file:
  
  res.frame <- read.table("monod.out.tmp",header=F)
  
  # return results:
  
  res <- res.frame[,2]
  names(res) <- res.frame[,1]   
  return(res)   
}


################################################################################


model.growth <- function(par,L)
{
  # model:
  # ------
  #
  #   growth of microorganisms on a substrate in a batch reactor:
  #
  #   dC_M         C_S
  #   ----  =  mu ----- C_M  -  b C_M
  #    dt         K+C_S
  #
  #   dC_S      mu   C_S
  #   ----  = - --  ----- C_M
  #    dt       Y   K+C_S
  #
  # state variables:
  # ----------------
  #
  #   C_M      concentration of microorganisms in the reactor
  #   C_S      concentration of substrate in the reactor
  #
  # parameters:
  # -----------
  #
  #   mu       maximum growth rate of microorganisms
  #   K        half-concentration of growth rate with respect to substrate
  #   b        rate of death and respiration processes of microorganisms
  #   Y        yield of growth process
  #   C_M_ini  initial concentration of microorganisms
  #   C_S_ini  initial concentration of substrate
  #
  # layout:
  # -------
  #
  #   L        observation layout containing variable names and evaluation time 
  #            points for the model
  #
  #            L should be a data frame with two columns. The first column 
  #            contains the names of the variables to be evaluated, and the 
  #            second column contains the time points (or other evaluation 
  #            steps) for the model. See help file for function "checkLayout" 
  #            in "bayesBias" package for more details. For example:
  #
  #              Monod model:           Growth model:
  #                  var  |   C            var  |   t
  #                ---------------       ---------------
  #                   r   |  0.1           C_M  |   1
  #                   r   |  0.2           C_M  |   2
  #                   r   |  0.3           C_M  |   3
  #                  ...  |  ...           ...  |  ...
  #                                        C_S  |   1
  #                                        C_S  |   2
  #                                        C_S  |   3
  #                                        ...  |  ...
  
  # check input:
  # ------------
  
  if ( length(par) < 6 )
  {
    stop("error in model.growth: insufficient number of parameters provided")
  }
  
  # process layout:
  # ---------------
  
  n <- nrow(L)
  L.names <- paste(L[,1],L[,2],sep="_")
  var <- unique(L[,1])
  if ( length(var) != 2 | is.na(match("C_M",var)) | is.na(match("C_S",var)) )
  {
    stop(paste(
      paste("*** error in model.growth:",
            "can only calculate results for variables \"C_M\" and \"C_S\";",
            "layout requests:"),
      paste(var,collapse=",")))
  }
  ind.M <- L[,1] == "C_M"
  ind.S <- L[,1] == "C_S"
  t <- sort(L[ind.M,2])
  if ( sum((sort(L[ind.S,2])-t)^2) != 0 )
  {
    stop(paste("*** error in model.growth:",
               "evalution times must be the same for C_M and C_S"))
  }
  
  # calculate results:
  # ------------------
  
  C_M_ini <- par["C_M_ini"]; if (is.na(C_M_ini)) C_M_ini <- par[5]  ## changed LW
  C_S_ini <- par["C_S_ini"]; if (is.na(C_S_ini)) C_S_ini <- par[6]  ## changed LW
  
  C_ini <- c(C_M=as.numeric(C_M_ini),C_S=as.numeric(C_S_ini))
  
  par   <- unlist(par)  ## changed LW
  
  res_ode <- ode(y=C_ini, times=t, func=rhs.growth, parms=par, method="rk4")
    ## changed LW
    ## using solver "ode" from "deSolve" package (method "rk4" is fast)

  # return results:
  # ---------------
  
  res_ode_df <- as.data.frame(res_ode)  ## changed LW (deSolve "ode" returns matrix)
  res <- as.numeric( c(res_ode_df$C_M, res_ode_df$C_S) )  ## changed LW
  names(res) <- c(paste("C_M_",t,sep=""),paste("C_S_",t,sep=""))
  
  # rearrange in required order:
  res.ord <- res[L.names]
  names(res.ord) <- L.names
  
  return(res.ord)
}


################################################################################


rhs.growth <- function(t,C,par)  ## changed LW (deSolve requires arguments in this order)
{
  # check input:
  # ------------
  
  if ( length(C) != 2 )
  {
    stop("error in rhs.growth: two state variables required")
  }
  if ( length(par) < 6 )
  {
    stop("error in rhs.growth: insufficient number of parameters provided")
  }
  
  # assign state variable and parameter values:
  # -------------------------------------------
  
  ind_M <- match("C_M",names(C)); if (is.na(ind_M)) ind_M <- 1
  ind_S <- match("C_S",names(C)); if (is.na(ind_S)) ind_S <- 2
  C_M   <- C[ind_M]
  C_S   <- C[ind_S]
  mu    <- par["mu"]; if (is.na(mu)) mu <- par[1]
  K     <- par["K"];  if (is.na(K))  K  <- par[2]
  b     <- par["b"];  if (is.na(b))  b  <- par[3]
  Y     <- par["Y"];  if (is.na(Y))  Y  <- par[4]
  
  # calculate results:
  # ------------------
  
  r_M <-         mu*C_S/(K+C_S)*C_M - b*C_M
  r_S <- - 1/Y * mu*C_S/(K+C_S)*C_M
  
  # return result:
  # --------------
  
  res <- numeric(0)
  res[ind_M] <- r_M
  res[ind_S] <- r_S
  return(list(res))  ## changed LW (deSolve "ode" requires as a list)
}


################################################################################