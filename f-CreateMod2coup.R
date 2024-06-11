#######################################################################################################
# Create POMP for SEIR with coupling
# See https://www.omnicalculator.com/physics/relative-humidity for details 
#######################################################################################################
CreateMod <- function(covars, lin_bool_val = T) {
  
  # Args: 
  # covars: table of covariates: week_no (starting at 0), Te and Td (in degrees celsius) (data frame)
  # lin_bool_val: should the transmission term be linearized? (boolean)
  # Returns: POMP model 
  
  # Check that covariates are correctly formatted
  # stopifnot(all(c("week_no", "Te", "Td", "RH_pred") %in% colnames(covars)))
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
  }"
    
  # Csnippet for deterministic skeleton
  det_process_model <- Csnippet("

    // Calculate transmission rate
    double beta11, beta12, beta21, beta22; 
    beta11 = beta22 = R0 * (mu + 1.0) / (1.0 + p); 
    beta12 = beta21 = beta11 * p; 

    // Climate-induced seasonality
    double RH_pred1 = pred_RH(Te1, Td1);
    double beta_seas1 = e_Te * (Te1 / Te_mean1 - 1.0) + e_RH * (RH_pred1 / RH_mean1 - 1.0); 
    beta_seas1 = exp(beta_seas1);
    
    double RH_pred2 = pred_RH(Te2, Td2);
    double beta_seas2 = e_Te * (Te2 / Te_mean2 - 1.0) + e_RH * (RH_pred2 / RH_mean2 - 1.0); 
    beta_seas2 = exp(beta_seas2);
    
    beta11 *= beta_seas1; 
    beta22 *= beta_seas2; 

    // Force of infection
    double lambda_t1 = beta11 * I1 / N1 + beta12 * I2 / N2; // Force of infection patch 1
    double lambda_t2 = beta21 * I1 / N1 + beta22 * I2 / N2; // Force of infection patch 2
    double p_SI1 = lin_bool ? lambda_t1 : (1.0 - exp(-lambda_t1)); // Proportion infected patch 1 during time step
    double p_SI2 = lin_bool ? lambda_t2 : (1.0 - exp(-lambda_t2)); // Proportion infected pathc 2 during time step

    // Difference equations
    DS1 = S1 +  mu * N1 + (1.0 - eps) * I1 + alpha * R1 - (p_SI1 + mu) * S1;
    DI1 = p_SI1 * S1 - mu * I1;
    DR1 = R1 + eps * I1 - (alpha + mu) * R1;
    DCC1 = p_SI1 * S1;

    DS2 = S2 +  mu * N2 + (1.0 - eps) * I2 + alpha * R2 - (p_SI2 + mu) * S2;
    DI2 = p_SI2 * S2 - mu * I2;
    DR2 = R2 + eps * I2 - (alpha + mu) * R2;
    DCC2 = p_SI2 * S2;
  ")
  
  # Csnippet for stochastic model
  stoch_process_model <- Csnippet("
  
    double S1_new, I1_new, R1_new, CC1_new, S2_new, I2_new, R2_new, CC2_new;
    
    // Calculate transmission rate
    double beta11, beta12, beta21, beta22; 
    beta11 = beta22 = R0 * (mu + 1.0) / (1.0 + p); 
    beta12 = beta21 = beta11 * p; 
    
    // Climate-induced seasonality
    double RH_pred1 = pred_RH(Te1, Td1);
    double beta_seas1 = e_Te * (Te1 / Te_mean1 - 1.0) + e_RH * (RH_pred1 / RH_mean1 - 1.0); 
    beta_seas1 = exp(beta_seas1);
    
    double RH_pred2 = pred_RH(Te2, Td2);
    double beta_seas2 = e_Te * (Te2 / Te_mean2 - 1.0) + e_RH * (RH_pred2 / RH_mean2 - 1.0); 
    beta_seas2 = exp(beta_seas2);
    
    beta11 *= beta_seas1; 
    beta22 *= beta_seas2; 
    
    // Force of infection
    double lambda_t1 = (beta11 * I1 / N1) + (beta12 * I2 / N2); // Force of infection patch 1
    double lambda_t2 = (beta21 * I1 / N1) + (beta22 * I2 / N2); // Force of infection patch 2
    
    // Environmental stochasticity
    double dW1 = rgammawn(sigma_beta, 1.0); 
    double dW2 = rgammawn(sigma_beta, 1.0);
    lambda_t1 *= dW1; 
    lambda_t2 *= dW2; 
    double p_SI1 = lin_bool ? lambda_t1 : (1.0 - exp(-lambda_t1)); // Proportion infected patch 1 during time step
    double p_SI2 = lin_bool ? lambda_t2 : (1.0 - exp(-lambda_t2)); // Proportion infected patch 2 during time step
    
    // Difference equations
    S1_new = S1 +  mu * N1 + (1.0 - eps) * I1 + alpha * R1 - (p_SI1 + mu) * S1; 
    I1_new = p_SI1 * S1 - mu * I1; 
    R1_new = R1 + eps * I1 - (alpha + mu) * R1; 
    CC1_new = p_SI1 * S1; 
    
    S2_new = S2 +  mu * N2 + (1.0 - eps) * I2 + alpha * R2 - (p_SI2 + mu) * S2; 
    I2_new = p_SI2 * S2 - mu * I2; 
    R2_new = R2 + eps * I2 - (alpha + mu) * R2; 
    CC2_new = p_SI2 * S2; 
    
    // Update
    S1 = S1_new; 
    I1 = I1_new; 
    R1 = R1_new; 
    CC1 = CC1_new; 
    
    S2 = S2_new; 
    I2 = I2_new; 
    R2 = R2_new; 
    CC2 = CC2_new; 
  ")
  
  # Csnippet for generating observations
  robs_model <- Csnippet("
    CC_obs1 = rnbinom_mu(1.0 / rho_k, rho_mean * CC1);
    CC_obs2 = rnbinom_mu(1.0 / rho_k, rho_mean * CC2);
  ")
  
  # Csnippet for evaluating likelihood of observations
  dobs_model <- Csnippet("
    double logL = ISNA(CC_obs1) ? 0.0 : dnbinom_mu(nearbyint(CC_obs1), 1.0 / rho_k, rho_mean * CC1, 1); // First component, on log-scale
    logL += ISNA(CC_obs2) ? 0.0 : dnbinom_mu(nearbyint(CC_obs2), 1.0 / rho_k, rho_mean * CC2, 1); // Second component, on log-scale
    lik = (give_log) ? logL : exp(logL);
  ")
  
  # Csnippet initial conditions 
  init_mod <- Csnippet("
    double S_star1 = N1 / R0; 
    S1 = S_star1; 
    I1 = ((alpha + mu) / (alpha + mu + eps) * (N1 - S_star1));
    R1 = eps / (alpha + mu + eps) * (N1 - S_star1);
    CC1 = 0; 
    
    double S_star2 = N2 / R0; 
    S2 = S_star2; 
    I2 = ((alpha + mu) / (alpha + mu + eps) * (N2 - S_star2));
    R2 = eps / (alpha + mu + eps) * (N2 - S_star2);
    CC2 = 0; 
    
    //S1 = N1-1; 
    //I1 = 1;
    //R1 = 0;
    //CC1 = 0; 
    
    //S2 = N2; 
    //I2 = 0;
    //R2 = 0;
    //CC2 = 0; 
  ")
  
  # Create model   
  mod <- pomp(
    data = data.frame(week = 1:(max(covars$week_no) - 1), 
                      CC_obs1 = NA,
                      CC_obs2 = NA), 
    times = "week", 
    t0 = 0, 
    obsnames = c("CC_obs1", "CC_obs2"), 
    covar = covariate_table(select(covars, c("week_no", "Te1", "Te2", "Td1", "Td2")), 
                            times = "week_no", order = "constant"), 
    statenames = c("S1", "I1", "R1", "CC1", "S2", "I2", "R2", "CC2"), 
    skeleton = pomp::map(f = det_process_model, delta.t = 1), 
    rprocess =  discrete_time(step.fun = stoch_process_model, delta.t = 1),
    rmeasure = robs_model, 
    dmeasure = dobs_model,
    rinit = init_mod,
    globals = RH_Cfun,
    params = c("mu" = 1 / 80 / 52, # Birth rate (per week)
               "p" = 0.01, # Coupling 
               "N1" = 5e6, # Total population size
               "N2" = 5e6, # Total population size
               "R0" = 2.5, # Reproduction no
               "sigma_beta" = 0.02, # SD of environmental noise
               "e_Te" = -0.2, # Effect of temperature
               "Te_mean1" = mean(covars$Te1, na.rm = T), # Temporal mean of Te 
               "Te_mean2" = mean(covars$Te2, na.rm = T), # Temporal mean of Te 
               "e_RH" = -0.2, # Effect of RH
               "RH_mean1" = mean(covars$RH_pred1, na.rm = T), # Temporal mean of RH
               "RH_mean2" = mean(covars$RH_pred2, na.rm = T), # Temporal mean of RH
               "eps" = 1, # Fraction of infections conferring immunity
               "alpha" = 1 / (2 * 52), # Waning rate 
               "rho_mean" = 0.1, # Average reporting probability
               "rho_k" = 0.04, # Reporting over-dispersion
               "lin_bool" = as.integer(lin_bool_val)), 
    paramnames = c("mu", "p", "N1", "N2", "R0", "sigma_beta", "e_Te", 
                   "Te_mean1", "Te_mean2", "e_RH", "RH_mean1", "RH_mean2", "eps", "alpha", 
                   "rho_mean", "rho_k", "lin_bool"),
    partrans = parameter_trans(log = c("R0", "alpha", "rho_k"), 
                               logit = c("rho_mean", "p"))
  )
  
  return(mod)
}

#######################################################################################################
# END
#######################################################################################################