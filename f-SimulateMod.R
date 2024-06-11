####################################################################################################
# Obtain simulations from the models in different locations 
####################################################################################################

fun_SimulateMod <- function(loc = "Cartagena",
                            parms = parms,
                            chg_parms = chg_parms,
                            nsim = 1) {
  # Args: 
  # loc: name of weather station for the target location (string)
  # pars: base parameters in the model (named vector)
  # chg_parms: parameters to change in the simulations (eps, R0, alpha with index) (data frame)
  # nsim: number of simulations (numeric)
  # Return: 
  # Table of simulated data
  
  # Create pomp model for a given location 
  PompMod <- CreateMod(covars = covars_l[[loc]], lin_bool_val = T)
  
  # Set parameters
  pomp::coef(PompMod, names(parms)) <- unname(parms)
  # Sets of changed parameters 
  p_mat <- parmat(params = pomp::coef(PompMod), nrep = nrow(chg_parms))
  p_mat["R0", ] <- chg_parms$R0
  p_mat["eps", ] <- chg_parms$eps
  p_mat["alpha", ] <- chg_parms$alpha
  
  # Run simulations
  sim_raw <- simulate(object = PompMod, nsim = nsim, params = p_mat, 
                      seed = 2186L, format = "data.frame")
  
  # If there is more than one simulation, but just one parameter set:
  if(nsim > 1 & length(chg_parms[[1]]) == 1) {sim_raw <- mutate(sim_raw, .id = paste0("1_", .id))}
  
  sim <- sim_raw %>% 
    separate(.id, into = c(".id", "sim_id"), sep = "_") %>%
    mutate(.id = as.integer(.id)) %>% 
    left_join(y = chg_parms) %>% 
    mutate(N_sim = S + I + R) %>% 
    rename(week_no = week) %>% 
    left_join(y = select(clim_dat_l[[loc]], week_no, matches("Te_norm|RH_pred_norm"), beta_seas)) 
  
  # If there is just one simulation: 
  if(nsim == 1) {sim <- select(sim, -sim_id)}
  
  # Data in long format
  sim_long <- sim %>% 
    pivot_longer(cols = c(S, I, R, N_sim, CC, CC_obs), names_to = "state_var", values_to = "value") %>% 
    mutate(state_var = factor(state_var, levels = c("S", "I", "R", "N_sim", "CC", "CC_obs")), 
           N = unname(pomp::coef(PompMod, "N")),
           loc = loc) 
  
  return(sim_long)
}
