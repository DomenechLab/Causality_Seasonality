#######################################################################################################
# Run simulations for vignette on quasi experiments
# Key point: different climates have different correlation between Te and RH, 
# making it more or less difficult to identify the effect of individual climatic variables
# In particular, all else equal tropical climates offer a "quasi-experiment" to identify the effect of RH
# Illustration: compare estimation performance in Rostock, Germany and Bogota, Colombia 
# In both locations, the individual variability of RH is comparable (CV=0.07 in Colombia, CV=0.08 in Rostock), 
# but the correlation with temperature is much higher in Rostock
# Hence, we expect easier identification of the effect of RH in Bogota and of Te in Rostock
#######################################################################################################

# Load packages -----------------------------------------------------------
rm(list = ls())
source("f-Pred_RH.R")
Pred_RH <- Vectorize(FUN = Pred_RH)
source("f-CreateMod.R")
source("f-CreateClimData.R")

source("s-base_packages.R")
library("pomp")
library("subplex")
library("corrplot")
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))
par(bty = "l", las = 1, lwd = 2)
save_plot <- F # Should all the plots be saved as a pdf? 

# Set model parameters ----------------------------------------------------
parms <- c("mu" = 1 / 80 / 52, # Birth rate 
           "N" = 5e6, # Total population size
           "R0" = 1.25, # Basic reproduction no
           "e_Te" = -0.2, # Effect of Te on transmission 
           "e_RH" = -0.2, # Effect of RH on transmission
           "eps" = 1, # Fraction of infections conferring sterilizing immunity (1: SIR, 0: SIS)
           "alpha" = 1 / (2 * 52), # Rate of waning immunity
           "rho_mean" = 0.1, # Average reporting probability 
           "rho_k" = 0.04) # Reporting over-dispersion

# Load climatic data in a given location------------------------------------------------------
# Bogota (Colombia): weather station "SKBO", rho(Te, RH) = -0.08, CV(RH) = 0.07
# Rostock (Germany): Rostock, rho(Te, RH) = -0.69, CV(RH) = 0.085

e_Te <- parms["e_Te"]
e_RH <- parms["e_RH"]
loc_nm <- "Rostock"
if(save_plot) pdf(file = sprintf("_saved/vignette-quasi-experiments-%s.pdf", loc_nm), width = 8, height = 8)

clim_dat <- CreateClimData(loc_nm = loc_nm, n_years = 10)

# Calculate seasonal term of transmission rate
clim_dat <- clim_dat %>% 
  mutate(beta_seas = exp(e_Te * (Te_norm - 1)) * exp(e_RH * (RH_pred_norm - 1)))

# Data in long format
clim_dat_long <- clim_dat %>% 
  pivot_longer(cols = Te:RH_pred_norm, names_to = "var", values_to = "value")

# Plot time series of climatic variables
pl <- ggplot(data = clim_dat_long %>% filter(var %in% c("Te", "Td", "RH", "RH_pred")), 
             mapping = aes(x = week_date, y = value)) + 
  geom_line() + 
  facet_wrap(~ var, scales = "free_y", ncol = 2) + 
  labs(x = "Time (weeks)", y = "Value", title = "Climatic data, time plot")
print(pl)

# Plot measured and predicted RH
pl <- ggplot(data = clim_dat, mapping = aes(x = RH, y = RH_pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  labs(title = "Agreement between RH and RH_pred")
print(pl)

# Plot time series of measured and predicted RH
pl <- ggplot(data = clim_dat, mapping = aes(x = week_date, y = RH)) + 
  geom_line() + 
  geom_line(mapping = aes(x = week_date, y = RH_pred), color = "red") + 
  labs(title = "Agreement between RH and RH_pred (red curve)")
print(pl)

# Season plots of all climatic variables
pl <- ggplot(data = clim_dat_long %>% filter(var %in% c("Te_norm", "Td_norm", "RH_pred_norm")), 
             mapping = aes(x = isoweek(week_date), y = value, 
                           #color = factor(year(week_date)), 
                           group = factor(year(week_date))
             )) + 
  geom_line(color = "grey") + 
  geom_smooth(mapping = aes(x = isoweek(week_date), y = value, group = NULL, color = NULL)) + 
  facet_wrap(~ var, scales = "fixed", ncol = 2) + 
  labs(x = "Week no", y = "Rescaled value", title = "Climatic data, season plot")
print(pl)

# Plot seasonal forcing of transmission rate
pl <- ggplot(data = clim_dat, mapping = aes(x = week_date, y = beta_seas)) + 
  geom_line() + 
  labs(x = "Time (weeks)", title = "Seasonal forcing in transmission rate")
print(pl)

# Add lagged variables for Te_norm and RH_norm
clim_dat <- clim_dat %>% 
  mutate(Te_norm_lag = lag(x = Te_norm, n = 1L, order_by = week_no), 
         RH_pred_norm_lag = lag(x = RH_pred_norm, n = 1L, order_by = week_no))

# Prepare covariate table
covars <- clim_dat %>% 
  select(week_date, week_no, Te, Td, RH_pred) %>% 
  arrange(week_date)

# Create pomp model -------------------------------------------------------
PompMod <- CreateMod(covars_df = covars, lin_bool_val = T)

# Set parameters
pomp::coef(PompMod, names(parms)) <- unname(parms)
base_pars <- pomp::coef(PompMod)
base_pars_trans <- pomp::coef(PompMod, transform = T)

# Run simulations ---------------------------------------------------------
pomp::coef(PompMod, names(base_pars)) <- unname(base_pars)
rho_mean_val <- unname(coef(PompMod, "rho_mean")) # Mean reporting probability
rho_k_val <- unname(coef(PompMod, "rho_k")) # Reporting over dispersion
N_val <- unname(coef(PompMod, "N"))


# Simulations
sim <- trajectory(object = PompMod, 
                  format = "data.frame") %>% 
  mutate(N_sim = S + I + R)

# Generate observations
sim$CC_obs <- freeze(expr = rnbinom(n = length(sim$CC), 
                                    mu = rho_mean_val * sim$CC, 
                                    size = 1 / rho_k_val), 
                     seed = 2186L)

# Pivot in long format
sim_long <- sim %>% 
  pivot_longer(cols = -c(".id", "week"), names_to = "state_var", values_to = "value") %>% 
  mutate(state_var = factor(state_var, levels = c("S", "I", "R", "N_sim", "CC", "CC_obs")), 
         N = N_val)

# Plot simulations
svars_plot <- c("S", "N_sim", "CC", "CC_obs")

pl <- ggplot(data = sim_long, 
             mapping = aes(x = week / 52, y = value)) + 
  geom_line() + 
  #scale_y_sqrt() +
  facet_wrap(~ state_var, scales = "free_y", ncol = 2, dir = "v") + 
  labs(x = "Year", y = "Proportion", title = "All model variables, time plot")
print(pl)

# Run estimation for climatic parameters -------------------------------------------------------------

n_rep <- 10 # No of replicate datasets

# Generate replicate datasets 
CC_obs_list <- freeze(seed = 2186L, 
                      expr = replicate(n = n_rep, 
                                       expr = rnbinom(n = length(sim$CC), 
                                                      mu = rho_mean_val * sim$CC, 
                                                      size = 1 / rho_k_val)))
fits <- vector(mode = "list", length = n_rep)

# Plot 
matplot(sim$week, CC_obs_list, type = "l", lty = 1, 
        xlab = "Week", ylab = "No of cases", 
        main = sprintf("%d generated datasets", n_rep))

# Parameters to estimate
all_pars_nm <- c("R0", "alpha", "rho_mean", "rho_k", "e_Te", "e_RH")

fits <- bake(file = sprintf("_saved/vignette-quasi-experiments-%s.rds", loc_nm), 
             expr = {
               for(s in seq_len(n_rep)) {
                 
                 print(s)
                 
                 # Assign data to pomp object
                 PompMod@data[1, ] <- CC_obs_list[, s]
                 
                 # Create objective function
                 LL_fun_cur <- traj_objfun(data = PompMod, est = all_pars_nm)
                 
                 # Run optimization
                 fits[[s]] <- try(
                   subplex(fn = LL_fun_cur, 
                           par = unname(base_pars_trans[all_pars_nm]), 
                           control = list(maxit = 1e6), 
                           hessian = T)
                 )
               }
               fits
             })

# Plot results ------------------------------------------------------------
base_pars_df <- data.frame(par = all_pars_nm, 
                           value = unname(base_pars_trans[all_pars_nm]))

# Remove estimations that failed,, if any
fits_class <- sapply(fits, class)
if(any(fits_class == "try-error")) fits <- fits[-which(fits_class == "try-error")]

# Extracts information on model fit
fits_info <- data.frame(sim_no = seq_along(fits), 
                        ll = -map_dbl(.x = fits, .f = "value"), 
                        conv = map_int(.x = fits, .f = "convergence"))

# Remove fits that didn't converge, if any
if(any(fits_info$conv < 0)) {
  fits <- fits[-which(fits_info$conv < 0)]
}

# Extract parameter estimates
pars_mle <- sapply(fits, getElement, "par")
rownames(pars_mle) <- all_pars_nm
pars_mle <- pars_mle %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(sim_no = 1:nrow(.)) %>% 
  select(sim_no, everything()) %>% 
  pivot_longer(cols = -sim_no, names_to = "par", values_to = "mle")

# Extract parameter SEs
pars_SE <- map_df(.x = fits, .f = function(x) {
  I <- solve(x$hessian)
  out <- sqrt(diag(I))
  names(out) <- all_pars_nm
  return(out)
})

pars_SE <- pars_SE %>% 
  mutate(sim_no = 1:nrow(.)) %>% 
  select(sim_no, everything()) %>% 
  pivot_longer(cols = -sim_no, names_to = "par", values_to = "se")

pars_all <- pars_mle %>% 
  left_join(y = pars_SE) %>% 
  left_join(y = base_pars_df %>% rename("true" = "value"))

# Calculate mean absolute bias 
est_perf <- pars_all %>% 
  group_by(par) %>% 
  summarise(mean_bias = mean(abs(mle - true)), 
            mean_se = mean(se)) %>% 
  ungroup() %>% 
  pivot_longer(cols = -par)

# Plot correlation matrix
R_mat <- fits[[1]]$hessian %>% 
  solve() %>% 
  cov2cor()
rownames(R_mat) <- colnames(R_mat) <- all_pars_nm

corrplot(corr = R_mat, method = "number", type = "lower")

# Plot confidence intervals for all parameters
pl <- ggplot(data = pars_all, 
             mapping = aes(x = mle, y = sim_no, xmin = mle - 2 * se, xmax = mle + 2 * se)) + 
  geom_pointrange() + 
  geom_vline(data = base_pars_df, mapping = aes(xintercept = value), linetype = "dashed", color = "red") + 
  facet_wrap(~ par, scales = "free_x") + 
  labs(x = "Estimate", y = "Simulation no", title = "Estimates, all parameters")
print(pl)

# Same plot for climatic parameters only
pl <- ggplot(data = pars_all %>% filter(par %in% c("e_Te", "e_RH")), 
             mapping = aes(x = mle, y = sim_no, xmin = mle - 2 * se, xmax = mle + 2 * se)) + 
  geom_pointrange() + 
  geom_vline(data = base_pars_df %>% filter(par %in% c("e_Te", "e_RH")), 
             mapping = aes(xintercept = value), linetype = "dashed", color = "red") + 
  facet_wrap(~ par, scales = "fixed", ncol = 1) + 
  labs(x = "Estimate", y = "Simulation no", title = "Estimates, climatic parameters")
print(pl)

# Plot estimation performance
pl <- ggplot(data = est_perf, mapping = aes(x = par, y = value)) + 
  geom_col() + 
  geom_text(color = "blue", mapping = aes(x = par, y = value + 0.01, label = round(value, 3))) + 
  facet_wrap(~ name, scales = "free_y", ncol = 1) + 
  labs(x = "Parameter", y = "Value", title = "Estimation performance")
print(pl)

# End  statements ---------------------------------------------------------------------
if(save_plot) dev.off()

#######################################################################################################
# END
#######################################################################################################
