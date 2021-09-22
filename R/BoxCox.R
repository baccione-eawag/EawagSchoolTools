################################################################################
#                                                                              #
# Functions for R package "bayesBias" : Bayesian Inference with Consideration  #
# of Bias                                                                      #
#                                                                              #
# Peter Reichert <peter.reichert@eawag.ch>                                     #
#                                                                              #
################################################################################


# BoxCox
# ======

# purpose:
# Box-Cox transformation

# arguments:
# data:        data to be transformed
# lambda1:     value of parameter lambda1 (exponent) of Box-Cox transformation
# lambda2:     value of parameter lambda2 (offset) of Box-Cox transformation

# output:
# transformed data

BoxCox <- function(data,lambda1=1,lambda2=1)
{
  if ( lambda1 == 0 )
  {
    return(ifelse(data>-lambda2,log(data+lambda2),NA))
  }
  else
  {
    return(ifelse(data>=-lambda2,((data+lambda2)^lambda1 - 1)/lambda1,NA))
  }
}


################################################################################


# BoxCox.deriv
# ============

# purpose:
# calculate derivative of Box-Cox transformation

# arguments:
# data:        data to be transformed
# lambda1:     value of parameter lambda1 (exponent) of Box-Cox transformation
# lambda2:     value of parameter lambda2 (offset) of Box-Cox transformation

# output:
# derivative of Box-Cox transformation

BoxCox.deriv <- function(data,lambda1=1,lambda2=1)
{
  return(ifelse(data>-lambda2,(data+lambda2)^(lambda1 - 1),NA))
}


################################################################################


# BoxCox.inv
# ==========

# purpose:
# inverse Box-Cox transformation

# arguments:
# data:        data to be back-transformed
# lambda1:     value of parameter lambda1 (exponent) of Box-Cox transformation
# lambda2:     value of parameter lambda2 (offset) of Box-Cox transformation

# output:
# back-transformed data

BoxCox.inv <- function(data,lambda1=1,lambda2=1)
{
  if ( lambda1 == 0 )
  {
    return(exp(data)-lambda2)
  }
  else
  {
    return(ifelse(lambda1*data>-1,(lambda1*data+1)^(1/lambda1)-lambda2,
                  -lambda2))
  }
}


################################################################################