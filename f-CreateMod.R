#######################################################################################################
# Calculate relative humidity (RH) from temperature and dew point temperature
# See https://www.omnicalculator.com/physics/relative-humidity for details 
#######################################################################################################
CreateMod <- function(covars_df, lin_bool_val = T) {
  
  # Args: 
  # covars_df: table of covariates: week_no (starting at 0), Te and Td (in degrees celsius) (data frame)
  # lin_bool_val: should the transmission term be linearized? (boolean)
  # Returns: POMP model 
  
  # Check that covariates are correctly formatted
  stopifnot(all(c("week_no", "Te", "Td", "RH_pred") %in% colnames(covars)))
  stopifnot(min(covars$week_no) == 0)
  covars <- arrange(covars, week_no)
  
  # C function for predicting RH from Te and Td
  RH_Cfun <- "
static double pred_RH(double Te, double Td) {
  // Args: 
  // Te: temperature (in degrees celsius)
  // Td: dew point temperature (in degrees celsius)
  // Returns: relative humidity (between 0 and 1)
  double beta_val = 17.625;
  double lambda_val = 243.04; 
  
  double out = beta_val * Td / (lambda_val + Td) - beta_val * Te / (lambda_val + Te);
  out = fmin2(1.0, exp(out));
  return out; 
}
"
  
  # Csnippet for deterministic skeleton
  det_process_model <- Csnippet("

  // Calculate transmission rate and force of infection
  double beta = R0 * (mu + 1.0); // Transmission rate
  double RH_pred = pred_RH(Te, Td);
  
  double beta_seas = e_Te * (Te / Te_mean - 1.0) + e_RH * (RH_pred / RH_mean - 1.0); 
  beta_seas = exp(beta_seas);
  double lambda_t = beta * I / N * beta_seas; // Force of infection
  double p_SI = lin_bool ? lambda_t : (1.0 - exp(-lambda_t)); // Proportion infected during time step
  
  // Difference equations
  DS = S +  mu * N + (1.0 - eps) * I + alpha * R - (p_SI + mu) * S; 
  DI = p_SI * S - mu * I; 
  DR = R + eps * I - (alpha + mu) * R; 
  DCC = p_SI * S; 
")
  
  # Csnippet for generating observations
  robs_model <- Csnippet("
  CC_obs = rnbinom_mu(1.0 / rho_k, rho_mean * CC);
")
  
  # Csnippet for evaluating likelihood of observations
  dobs_model <- Csnippet("
  lik = dnbinom_mu(nearbyint(CC_obs), 1.0 / rho_k, rho_mean * CC, give_log);
")
  
  # Csnippet for initializing state variables
  # All variables are initialized at the endemic equilibrium of the seasonally-unforced model
  init_mod <- Csnippet("
  double S_star = N / R0; 
  S = S_star; 
  I = (alpha + mu) / (alpha + mu + eps) * (N - S_star);
  R = eps / (alpha + mu + eps) * (N - S_star);
  CC = 0; 
")
  
  # Create pomp object
  mod <- pomp(data = data.frame(week = 1:(max(covars$week_no) - 1), CC_obs = NA), 
              times = "week", 
              t0 = 0, 
              obsnames = "CC_obs", 
              covar = covariate_table(select(covars, c("week_no", "Te", "Td")), 
                                      times = "week_no", order = "constant"), 
              statenames = c("S", "I", "R", "CC"), 
              skeleton = pomp::map(f = det_process_model, delta.t = 1),  
              rmeasure = robs_model, 
              dmeasure = dobs_model,
              rinit = init_mod,
              globals = RH_Cfun,
              params = c("mu" = 1 / 80 / 52, # Birth rate (per week)
                         "N" = 5e6, # Total population size
                         "R0" = 1.25, # Reproduction no
                         "e_Te" = 0, # Effect of temperature
                         "Te_mean" = mean(covars$Te, na.rm = T), # Temporal mean of Te 
                         "e_RH" = 0, # Effect of RH
                         "RH_mean" = mean(covars$RH_pred, na.rm = T), # Temporal mean of RH
                         "eps" = 1, # Fraction of infections conferring immunity
                         "alpha" = 0, # Waning rate 
                         "rho_mean" = 1, # Average reporting probability
                         "rho_k" = 0.04, # Reporting over-dispersion
                         "lin_bool" = as.integer(lin_bool_val)), 
              paramnames = c("mu", "N", "R0", "e_Te", "Te_mean", 
                             "e_RH", "RH_mean", "eps", "alpha", "rho_mean",
                             "rho_k", "lin_bool"), 
              partrans = parameter_trans(log = c("R0", "alpha", "rho_k"), 
                                         logit = c("rho_mean"))
  )
  
  return(mod)
}