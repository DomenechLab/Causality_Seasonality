#######################################################################################################
# Calculate relative humidity (RH) from temperature and dew point temperature
# See https://www.omnicalculator.com/physics/relative-humidity for details 
#######################################################################################################

Pred_RH <- function(temp, dewPoint) {
  # Args: 
  # temp: temperature (in degrees celsius) (numeric)
  # dewPoint: dew point temperature (in degrees celsius) (numeric)
  # Returns: relative humidity (between 0 and 1) (numeric)
  
  beta <- 17.625
  lambda <- 243.04
  
  out <- beta * dewPoint / (lambda + dewPoint) - beta * temp / (lambda + temp)
  out <- min(1, exp(out))
  return(out)
}

####################################################################################################
# End
####################################################################################################