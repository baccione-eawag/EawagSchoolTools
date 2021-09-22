################################################################################
#                                                                              #
# Functions for R package "eawagSummerSchoolTools" : Miscellaneous Tools for   #
# Eawag Summer School in Environmental Systems Analysis                        #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# Gnls
# ====

# purpose:
# calculate the generalized least squares parameter estimates of a nonlinear model

# arguments:
# model.name:  name of the function representing the model
#              (note that this model requires the parameters to be specified
#              explicitly in the function headers)
# y.obs:       observed data corresponding to model output
# var.inv:     inverse variance-covariance matrix of the multivariate normal distribution
#              that characterizes the error model
# par:         named parameter vector with initial values; 
#              passed to the model as separate arguments
# x:           list of named inputs passed to the model
# ...          further optional arguments passed to nls

# output:
# list of:
# model.name:  name of the function representing the model
# par:         named parameter vector with parameter estimates 
# y.obs:       observed data corresponding to model output
# y.det:       deterministic model results corresponding to the estimated parameters
# resid:       residuals
# var.inv:     inverse variance-covariance matrix of the multivariate normal distribution
#              that characterizes the error model
# x:           input object passed as first argument to the model
# res.nls:     results from application of the nls function

Gnls <- function(model.name,y.obs,var.inv,par,x=list(),...)
{
  A <- chol(var.inv)
  add.args <- ""
  if ( names(x) > 0 )
  {
    add.args <- paste(",",paste(names(x),collapse=","),sep="")
  }
  model.call <- paste(model.name,"(",paste(names(par),collapse=","),
                      add.args,")",sep="")
  data <- x
  data$y.obs.trans <- A%*%y.obs
  res.nls <- nls(as.formula(paste("y.obs.trans ~ A%*%",model.call)),
                 data=data,
                 start=par,...)
  
  coef  <- coef(res.nls)
  args <- as.list(coef)
  if ( length(x) > 0 ) args <- c(args,x)
  y.det <- do.call(model.name,args=args)  # reevaluate to avoid need for
  resid <- y.obs - y.det                  # back transformation
  
  return(list(model.name = model.name,
              par        = coef,
              y.obs      = y.obs,
              y.det      = y.det,
              resid      = resid,
              var.inv    = var.inv,
              x          = x,
              res.nls    = res.nls))
}


################################################################################


# Gnls.diag
# =========

# purpose:
# calculate regression diagnostics for generalized least squares

# arguments:
# res.gnls:    output from Gnls
# par.inc:     increments of parameters used to approximate the derivatives

# output:
# list of:
# par.inc:     parameter increments used to calculate V
# V:           matrix of partial derivatives of model results with respect to
#              parameters at the estimated parameter values
# sd.par:      standard errors of the parameters assuming known error variances
# var.par:     estimated variance-covariance matrix of the parameters 
#              assuming known error variances
# corr.par:    estimated correlation matrix of the parameters
# sd.rel:      estimated correction factor in standard deviations of the error model
# sd.par.rel:  standard errors of the parameters when estimating a common factor
#              in error variances
# var.par.rel: estimated variance-covariance matrix of the parameters 
#              when estimating a common factor in error variances

Gnls.diag <- function(res.gnls,par.inc)
{
  # sensitivites:
  # -------------
  
  V <- sens.loc.explicitpars(par        = res.gnls$par,
                             model.name = res.gnls$model.name,
                             par.inc    = par.inc,
                             x          = res.gnls$x)
  
  # calculate variance-covariance matrix of the estimator for given error variance:
  # -------------------------------------------------------------------------------
  
  var.par  <- solve( t(V) %*% res.gnls$var %*% V )
  rownames(var.par) = names(res.gnls$par)
  colnames(var.par) = names(res.gnls$par)
  
  sd.par   <- sqrt(diag(var.par))
  names(sd.par) <- names(res.gnls$par)
  
  corr.par <- (1/sd.par) %*% t(1/sd.par) * var.par
  rownames(corr.par) = names(res.gnls$par)
  colnames(corr.par) = names(res.gnls$par)
  
  # variance-covariance matrix of the estimator for estimated error variance:
  # -------------------------------------------------------------------------
  
  sd.rel <- sqrt(as.numeric(
    t(res.gnls$y.obs-res.gnls$y.det) %*% res.gnls$var.inv %*% (res.gnls$y.obs-res.gnls$y.det) /
      (length(res.gnls$y.obs)-length(res.gnls$par))
  ))
  var.par.rel <- var.par*sd.rel^2
  sd.par.rel  <- sd.par*sd.rel
  
  return(c(res.gnls,
           list(par.inc     = par.inc,
                V           = V,
                sd.par      = sd.par,
                var.par     = var.par,
                corr.par    = corr.par,
                sd.rel      = sd.rel,
                sd.par.rel  = sd.par.rel,
                var.par.rel = var.par.rel
           )))
}


################################################################################


# Gnls.test
# =========

# purpose:
# calculate test statistics for generalized least squares

# arguments:
# res.gnls:    output from Gnls

# output:
# list of:
# V:           matrix of partial derivatives of model results with respect to
#              parameters at the estimated parameter values
# sd.par:      standard errors of the parameters assuming known error variances
# var.par:     estimated variance-covariance matrix of the parameters 
#              assuming known error variances
# corr.par:    estimated correlation matrix of the parameters
# sd.rel:      estimated correction factor in standard deviations of the error model
# sd.par.rel:  standard errors of the parameters
#              when estimating a common factor in error variances
# var.par.rel: estimated variance-covariance matrix of the parameters 
#              when estimating a common factor in error variances

Gnls.test <- function(res.gnls,V,par,y.det)
{
  # chi2 test; to be compared with 
  # qchisq(1-alpha,length(y.obs)):
  # ------------------------------
  
  chi2 <- t(res.gnls$y.obs-y.det) %*% res.gnls$var.inv %*% (res.gnls$y.obs-y.det)
  
  # exact F test; to be compared with 
  # qf(1-alpha,length(par),length(y.obs)-length(par)):
  # --------------------------------------------------
  
  F.exact <- t(res.gnls$y.obs-y.det) %*% res.gnls$var.inv %*% V %*% 
    solve(t(V)%*%res.gnls$var.inv%*%V) %*% 
    t(V) %*% res.gnls$var.inv %*% (res.gnls$y.obs-y.det) /
    t(res.gnls$y.obs-y.det) %*% 
    ( res.gnls$var.inv -
        res.gnls$var.inv %*% V %*% 
        solve(t(V)%*%res.gnls$var.inv%*%V) %*% 
        t(V) %*% res.gnls$var.inv ) %*% 
    (res.gnls$y.obs-y.det) *
    (length(res.gnls$y.obs)-length(res.gnls$par))/length(res.gnls$par)
  
  # F test for linearized model; to be compared with 
  # qf(1-alpha,length(par),length(y.obs)-length(par)):
  # --------------------------------------------------
  
  F.lin <- t(par-res.gnls$par) %*% t(V) %*% res.gnls$var.inv %*% V  %*% (par-res.gnls$par) /
    ( t(res.gnls$y.obs-res.gnls$y.det) %*% res.gnls$var.inv %*% (res.gnls$y.obs-res.gnls$y.det) ) *
    (length(res.gnls$y.obs)-length(res.gnls$par))/length(res.gnls$par)
  
  return(list(chi2        = chi2,
              F.exact     = F.exact,
              F.lin       = F.lin))
}


################################################################################


# Gnls.predict
# ============

# purpose:
# calculate predictions based on generalized least squares regression results

# arguments:
# diag.gnls:   output from Gnls.diag
# newx:        new model input x at which confidence intervals are to be calculated
# level:       probability level of conficence intervals (default 0.95)
# var:         error variance-covariance matrix at new model input (only variances used)

# output:
# list of:
# x:           x values at which model results are calculated
# y.pred:      predicted model results
# confint:     confidence intervals of deterministic model results
# confint.rel: confidence intervals of deterministic model results
#              when estimating a common factor in error variances
# predint:     prediction intervals (only if error variance.covariance matrix is provided)
# predint.rel: prediction intervals (only if error variance.covariance matrix is provided)
#              when estimating a common factor in error variances

Gnls.predict <- function(diag.gnls,newx,level=0.95,var=NA)
{
  # prediction and sensitivites:
  # ----------------------------
  
  args <- as.list(diag.gnls$par)
  if ( length(newx) > 0 ) args <- c(args,newx)
  res.par <- do.call(diag.gnls$model.name,args=args)
  newV <- sens.loc.explicitpars(par        = diag.gnls$par,
                                model.name = diag.gnls$model.name,
                                par.inc    = diag.gnls$par.inc,
                                x          = newx)
  
  # variance-covariance matrix and standard deviaitons of model results:
  # --------------------------------------------------------------------
  
  var.y    <- newV %*% diag.gnls$var.par %*% t(newV)
  sd.y     <- sqrt(diag(var.y))
  sd.y.rel <- diag.gnls$sd.rel*sd.y
  
  # calculate confidence intervals:
  # -------------------------------
  
  df          <- length(diag.gnls$y.obs)-length(diag.gnls$par)
  confint     <- Confint(est=res.par,sd=sd.y,df=NA,level=level)
  confint.rel <- Confint(est=res.par,sd=sd.y.rel,df=df,level=level)
  
  res <- list(x           = newx,
              y.pred      = res.par,
              confint     = confint,
              confint.rel = confint.rel)
  
  if ( is.matrix(var) )
  {
    sd.y.err     <- sqrt(sd.y^2+diag(var))
    sd.y.rel.err <- sqrt(sd.y.rel^2+diag(var)*diag.gnls$sd.rel^2)
    predint     <- Confint(est=res.par,sd=sd.y.err,df=NA,level=level)
    predint.rel <- Confint(est=res.par,sd=sd.y.rel.err,df=df,level=level)
    res <- c(res,list(predint=predint,predint.rel=predint.rel))
  }
  
  return(res)
}


################################################################################