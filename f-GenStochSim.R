#######################################################################################################
# Generate stochastic simulation for pomp model
# Stochasticity is introduced via environmental noise in the transmission rate
# Everything else is assumed deterministic
#######################################################################################################

GenStochsim <- function(pomp_mod, beta_sigma = 0.05) {
  
  # Args: 
  # pomp_mod: POMP model object
  # beta_sigma: SD of environmental Gamma white noise in transmission
  # Returns: data frame with state variables
  
  # Extract pomp components
  pars <- coef(pomp_mod) # Model parameters
  times_mod <- pomp_mod@times # Simulation times
  times_mod <- c(0, times_mod) # Add time 0
  x0 <- rinit(pomp_mod)
  x0 <- setNames(object = as.numeric(x0), nm = rownames(x0)) # Initial conditions
  x0 <- c(x0, "CC_obs" = 0) # Add observation variable
  covars_df <- pomp_mod@covar@table %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(week = pomp_mod@covar@times) %>% 
    select(week, everything())
  rownames(covars_df) <- NULL
  
  # Create output data frame
  out <- matrix(data = NA, 
                nrow = length(times_mod), 
                ncol = length(x0) + 1, 
                dimnames = list(NULL, c("week", names(x0)))) %>% 
    as.data.frame()
  
  # Initialize
  out[1, names(x0)] <- unname(x0)
  out$week <- times_mod
  
  for(tt in 2:nrow(out)) {
    
    t_cur <- out$week[tt - 1]
    S_cur <- out$S[tt - 1]
    I_cur <- out$I[tt - 1]
    R_cur <- out$R[tt - 1]
    Te_cur <- covars_df$Te[covars_df$week == (tt - 1)]
    Td_cur <- covars_df$Td[covars_df$week == (tt - 1)]
    
    #Calculate transmission rate
    
    # Climate-induced seasonal forcing
    RH_pred = Pred_RH(temp = Te_cur, dewPoint = Td_cur)
    beta_seas <-  pars["e_Te"] * (Te_cur / pars["Te_mean"] - 1) + pars["e_RH"] * (RH_pred / pars["RH_mean"] - 1)
    beta_seas <-  exp(beta_seas) # Seasonal forcing
    b <-  pars["R0"] * (pars["mu"] + 1.0) * beta_seas # Transmission rate
    
    # Environmental stochasticity
    dW <-  rgammawn(n = 1, sigma = beta_sigma, dt = 1.0) # Mean: 1, SD: beta_sigma 
    b <- b * dW 
    
    # Calculate force of infection
    lambda_t <-  b * I_cur / pars["N"]  #Force of infection
    p_SI <-  ifelse(pars["lin_bool"], lambda_t, (1.0 - exp(-lambda_t))) # Proportion infected during time step
    
    # Difference equations
    S_new <-  S_cur +  pars["mu"] * pars["N"] + (1 - pars["eps"]) * I_cur + pars["alpha"] * R_cur - (p_SI + pars["mu"]) * S_cur 
    I_new <-  p_SI * S_cur - pars["mu"] * I_cur 
    R_new <-  R_cur + pars["eps"] * I_cur - (pars["alpha"] + pars["mu"]) * R_cur
    CC_new <-  p_SI * S_cur
    
    # Update
    out[tt, c("S", "I", "R", "CC")] <- c(S_new, I_new, R_new, CC_new)
  }
  
  # Generate observation model
  out$CC_obs <-  rnbinom(n = length(out$CC), mu = pars["rho_mean"] * out$CC, size = 1 / pars["rho_k"])
  
  # Remove time 0
  out <- filter(out, week > 0)
  
  return(out)
}