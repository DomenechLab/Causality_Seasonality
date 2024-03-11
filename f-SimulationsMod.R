#######################################################################################################
# Obtain simulations from the models in different locations 
#######################################################################################################

fun_SimulationsMod <- function(loc = "Cartagena",
                               all_parms = all_parms,
                               nsim = 1) {
  # Args: 
  # loc: name of weather station for the target location (string)
  # all_params: parameters to change in the simulations (eps, R0, alpha with index) (data frame)
  # nsim: number of simulations (numeric)
  # Return: 
  # Table of absolute and normalized values for weekly means for Te, Td, RH, and RH_pred (data frame)
  
  # Create pomp model for a given location 
  PompMod <- CreateMod(covars = covars_l[[loc]], lin_bool_val = T)
  
  # Set parameters
  pomp::coef(PompMod, names(parms)) <- unname(parms)
  base_pars <- coef(PompMod)
  
  # Run simulations
  pomp::coef(PompMod, names(base_pars)) <- unname(base_pars)
  rho_mean_val <- unname(pomp::coef(PompMod, "rho_mean"))
  #rho_k_val <- 0.04 # Reporting overdispersion
  N_val <- unname(pomp::coef(PompMod, "N"))
  
  p_mat <- parmat(params = pomp::coef(PompMod), 
                  nrep = nrow(all_parms))
  p_mat["eps", ] <- all_parms$eps
  p_mat["R0", ] <- all_parms$R0
  p_mat["alpha", ] <- all_parms$alpha
  
  # Run trajectory and generate observations 
  sim_raw <- simulate(object = PompMod, 
                      nsim = nsim,
                      params = p_mat, 
                      format = "data.frame")
  
  # If there is more than one simulation, but just one parameter set:
  if(nsim > 1 & length(all_parms[[1]]) == 1) {
    sim_raw <- mutate(sim_raw, .id = paste0("1_", .id))
    }
  
  sim <- sim_raw %>% 
    separate(.id, into = c(".id", "sim_id"), sep = "_") %>%
    mutate(.id = as.integer(.id)) %>% 
    left_join(y = all_parms) %>% 
    mutate(N_sim = S + I + R) %>% 
    #CC_obs = rnbinom(n = length(CC), mu = rho_mean_val * CC, size = 1 / rho_k_val)) %>% 
    rename(week_no = week) %>% 
    left_join(y = clim_dat_l[[loc]] %>% 
                select(week_no, matches("Te_norm|RH_pred_norm"), beta_seas)) %>% 
    select(.id, sim_id, eps, R0, alpha, week_no, S, I, R, N_sim, CC, CC_obs, everything())
  
  # If there is just one simulation: 
  if(nsim == 1) {sim <- select(sim, -sim_id)}
  
  # Data in long format
  sim_long <- sim %>% 
    pivot_longer(cols = S:CC_obs, names_to = "state_var", values_to = "value") %>% 
    mutate(state_var = factor(state_var, levels = c("S", "I", "R", "N_sim", "CC", "CC_obs")), 
           N = N_val,
           loc = loc) 
  
  return(sim_long)
}
