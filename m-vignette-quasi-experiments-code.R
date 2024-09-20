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
source("f-GenStochSim.R")

source("s-base_packages.R")
library("bbmle")
library("pomp")
library("corrplot")

theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))
par(bty = "l", las = 1, lwd = 2)
save_plot <- F # Should all the plots be saved as a pdf? 

# Set model parameters ----------------------------------------------------
parms <- c("mu" = 1 / 80 / 52, # Birth rate 
           "N" = 5e6, # Total population size
           "R0" = 2.5, # Basic reproduction no
           "sigma_beta" = 0.02, # SD of environmental noise (set to 0 for a deterministic process model)
           "e_Te" = -0.2, # Effect of Te on transmission 
           "e_RH" = -0.2, # Effect of RH on transmission
           "eps" = 1, # Fraction of infections conferring sterilizing immunity (1: SIR, 0: SIS)
           "alpha" = 1 / (1 * 52), # Rate of waning immunity
           "rho_mean" = 0.1, # Average reporting probability 
           "rho_k" = 0.04) # Reporting over-dispersion

# Load climatic data in a given location------------------------------------------------------
# Bogota (Colombia): weather station "SKBO", rho(Te, RH) = -0.08, CV(RH) = 0.07
# Rostock (Germany): Rostock, rho(Te, RH) = -0.69, CV(RH) = 0.085
# Pasto (Colombia): "SKPS", rho(Te, RH) = -91, CV(RH) = 0.019

e_Te <- parms["e_Te"]  
e_RH <- parms["e_RH"]
loc_nm <- "SKBO"
#if(save_plot) pdf(file = sprintf("_figures/vignette-quasi-experiments-%s.pdf", loc_nm), width = 8, height = 8)

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
  labs(x = "Time (weeks)", y = "Value", 
       title = sprintf("Climatic data, time plot (location: %s)", loc_nm))
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
PompMod <- CreateMod(covars = covars, lin_bool_val = T)

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
  labs(x = "Year", y = "Proportion", 
       title = "All model variables, time plot (deterministic model)", 
       subtitle = sprintf("R0=%.2f, 1/alpha=%.1f yr, rho_mean=%.1f, rho_k=%.2f\nN=%.1f M, eps=%.1f, e_Te=%.1f, e_RH=%.1f, s_beta=%.2f", 
                          parms["R0"], 1 / parms["alpha"] / 52, parms["rho_mean"], parms["rho_k"], parms["N"] / 1e6, parms["eps"], 
                          parms["e_Te"], parms["e_RH"], parms["sigma_beta"]))
print(pl)

# Generate observations, assuming only observation noise ---------------------------------------------------
n_rep <- 100 # No of replicate datasets

# Generate replicate datasets 
CC_obs_noiseObs <- freeze(seed = 2186L, 
                          expr = replicate(n = n_rep, 
                                           expr = rnbinom(n = length(sim$CC), 
                                                          mu = rho_mean_val * sim$CC, 
                                                          size = 1 / rho_k_val)))

# Plot 
matplot(sim$week, CC_obs_noiseObs, type = "l", lty = 1, 
        xlab = "Week", ylab = "No of cases", 
        main = sprintf("%d generated datasets, with observation noise only", n_rep), 
        col = "grey")
lines(sim$week, sim$CC * parms["rho_mean"], col = "red", lwd = 2)

# Generate observation, assuming observation AND process noise ----------------------------------------------------

sim_noiseAll <- simulate(object = PompMod, 
                         nsim = n_rep, 
                         seed = 2186L, 
                         format = "data.frame")

CC_obs_noiseAll <- sim_noiseAll %>% 
  rename("sim" = ".id") %>% 
  select(sim, week, CC_obs) %>% 
  pivot_wider(names_from = "sim", values_from = "CC_obs") %>% 
  arrange(week) %>% 
  select(-week) %>% 
  as.data.frame()

matplot(sim$week, CC_obs_noiseAll, type = "l", lty = 1, 
        xlab = "Week", ylab = "No of cases", 
        main = sprintf("%d generated datasets, with observation and process noise", n_rep), 
        col = "grey")
lines(sim$week, sim$CC * parms["rho_mean"], col = "red", lwd = 2)

# Run estimation for climatic parameters -------------------------------------------------------------

fits <- vector(mode = "list", length = n_rep)

# Parameters to estimate
all_pars_nm <- c("R0", "alpha", "rho_mean", "rho_k", "e_Te", "e_RH")

fits <- bake(file = paste0("_saved/_vignette-quasi-experiments/R0_", parms[["R0"]],"/vignette-quasi-experiments-", loc_nm, ".rds"), 
             expr = {
               for(s in seq_len(n_rep)) {
                 
                 print(s)
                 
                 # Assign data to pomp object
                 PompMod@data[1, ] <- as.numeric(CC_obs_noiseAll[, s])
                 
                 # Create objective function
                 LL_fun_cur <- traj_objfun(data = PompMod, est = all_pars_nm)
                 parnames(LL_fun_cur) <- all_pars_nm
                 
                 # Run optimization
                 fits[[s]] <- try(
                   mle2(minuslogl = LL_fun_cur, 
                        start = base_pars_trans[all_pars_nm], 
                        optimizer = "optim", 
                        method = "Nelder-Mead", 
                        control = list(maxit = 1e6))
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
                        counts = map_dbl(.x = fits, .f = ~ .x@details$counts[1]),
                        ll = -map_dbl(.x = fits, .f = ~ .x@details$value), 
                        conv = map_int(.x = fits, .f = ~ .x@details$convergence))

# Remove fits that didn't converge, if any
if(any(fits_info$conv < 0)) {
  fits <- fits[-which(fits_info$conv < 0)]
}

# Extract MLE and SE of estimates
pars_mle_se <- purrr::map(.x = fits, 
                          .f = ~ .x %>% summary() %>% coef() %>% as.data.frame() %>% rownames_to_column())

pars_mle_se <- pars_mle_se %>% 
  bind_rows(.id = "sim_no") %>% 
  rename("par" = "rowname", 
         "mle" = "Estimate", 
         "se" = "Std. Error") %>% 
  select(sim_no, par, mle, se) %>% 
  left_join(y = base_pars_df %>% rename("true" = "value")) %>% 
  mutate(sim_no = as.integer(sim_no)) %>% 
  left_join(y = fits_info)

# Replace NaN SE with NAs
if(any(is.nan(pars_mle_se$se))) {
  pars_mle_se$se[is.nan(pars_mle_se$se)] <- NA
}

# Calculate mean absolute bias 
est_perf <- pars_mle_se %>% 
  group_by(par) %>% 
  summarise(mean_bias = mean(abs(mle - true)), 
            mean_se = mean(se, na.rm = T)) %>% 
  ungroup() %>% 
  pivot_longer(cols = -par)

# Plots -------------------------------------------------------------------
# Plot correlation matrix
# R_mat <- fits[[1]]@vcov %>% cov2cor()
# 
# corrplot(corr = R_mat, method = "number", type = "lower")

# Plot confidence intervals for all parameters
pl <- ggplot(data = pars_mle_se, 
             mapping = aes(x = mle, y = sim_no, xmin = mle - 2 * se, xmax = mle + 2 * se, color = ll)) + 
  geom_pointrange() + 
  geom_vline(data = base_pars_df, mapping = aes(xintercept = value), linetype = "dashed", color = "red") + 
  facet_wrap(~ par, scales = "free_x") + 
  scale_color_viridis(option = "viridis") + 
  labs(x = "Estimate", y = "Simulation no", title = "Estimates, all parameters")
print(pl)

# Same plot for climatic parameters only
pl <- ggplot(data = pars_mle_se %>% filter(par %in% c("e_Te", "e_RH")), 
             mapping = aes(x = mle, y = sim_no, xmin = mle - 2 * se, xmax = mle + 2 * se, color = ll)) + 
  geom_pointrange() + 
  geom_vline(data = base_pars_df %>% filter(par %in% c("e_Te", "e_RH")), 
             mapping = aes(xintercept = value), linetype = "dashed", color = "red") + 
  facet_wrap(~ par, scales = "fixed", ncol = 1) + 
  scale_color_viridis(option = "viridis") + 
  labs(x = "Estimate", y = "Simulation no", title = "Estimates, climatic parameters")
print(pl)

# Plot estimation performance
pl <- ggplot(data = est_perf, mapping = aes(x = par, y = value)) + 
  geom_col() + 
  geom_text(color = "blue", mapping = aes(x = par, y = value + 0.01, label = round(value, 3))) + 
  facet_wrap(~ name, scales = "free_y", ncol = 1) + 
  labs(x = "Parameter", y = "Value", title = "Estimation performance")
print(pl)

# Save everything -------------------------------------------------------
to_save <- list(sim_det = sim, 
                sim_stoch = sim_noiseAll, 
                pars_est = pars_mle_se)

saveRDS(object = to_save, 
        file = paste0("_saved/_vignette-quasi-experiments/R0_", parms[["R0"]],"/vignette-quasi-experiments-", loc_nm, "-all.rds")) 

# End  statements ---------------------------------------------------------------------
if(save_plot) dev.off()

#######################################################################################################
# End
#######################################################################################################
