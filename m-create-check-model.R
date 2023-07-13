#######################################################################################################
# Run simulations to check the model
#######################################################################################################

# Load packages -----------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
source("f-Pred_RH.R")
Pred_RH <- Vectorize(FUN = Pred_RH)
source("f-CreateClimData.R")
source("f-CreateMod.R")
library(pomp)
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))

# Set model parameters ----------------------------------------------------
parms <- c("mu" = 1 / 80 / 52, # Birth rate 
           "N" = 5e6, # Total population size
           "R0" = 1.25, # Basic reproduction no
           "sigma_beta" = 0, # SD of environmental noise (set to 0 for a deterministic process model)
           "e_Te" = -0.2, # Effect of Te on transmission 
           "e_RH" = -0.2, # Effect of RH on transmission
           "eps" = 1, # Fraction of infections conferring sterilizing immunity (1: SIR, 0: SIS)
           "alpha" = 1 / (1 * 52), # Rate of waning immunity (per week)
           "rho_mean" = 0.1, # Average reporting probability 
           "rho_k" = 0.04) # Reporting over-dispersion

# Load climatic data in a given location------------------------------------------------------
# Bogota: weather station "SKBO"
# Madrid: weather station "LEVS"
clim_dat <- readRDS("_data/clim_data.rds")

e_Te <- parms["e_Te"]
e_RH <- parms["e_RH"]

clim_dat <- CreateClimData(loc_nm = "Rostock", n_years = 10)

# Calculate seasonal term of transmission rate
clim_dat <- clim_dat %>% 
  mutate(beta_seas = exp(e_Te * (Te_norm - 1)) * exp(e_RH * (RH_pred_norm - 1)))

clim_dat_long <- clim_dat %>% 
  pivot_longer(cols = Te:RH_pred_norm, names_to = "var", values_to = "value")

pl <- ggplot(data = clim_dat_long %>% filter(var %in% c("Te", "Td", "RH")), 
             mapping = aes(x = week_date, y = value)) + 
  geom_line() + 
  facet_wrap(~ var, scales = "free_y", ncol = 2)
print(pl)

# Plot measured and predicted RH
pl <- ggplot(data = clim_dat, mapping = aes(x = RH, y = RH_pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")
print(pl)

# Plot time series of measured and predicted RH
pl <- ggplot(data = clim_dat, mapping = aes(x = week_date, y = RH)) + 
  geom_line() + 
  geom_line(mapping = aes(x = week_date, y = RH_pred), color = "red")
print(pl)

# Season plots of all climatic variables
pl <- ggplot(data = clim_dat_long %>% filter(var %in% c("Te_norm", "Td_norm", "RH_norm")), 
             mapping = aes(x = isoweek(week_date), y = value, 
                           #color = factor(year(week_date)), 
                           group = factor(year(week_date))
             )) + 
  geom_line(color = "grey") + 
  geom_smooth(mapping = aes(x = isoweek(week_date), y = value, group = NULL, color = NULL)) + 
  facet_wrap(~ var, scales = "fixed", ncol = 2)
print(pl)

# Plot seasonal forcing of transmission rate
pl <- ggplot(data = clim_dat, mapping = aes(x = week_date, y = beta_seas)) + 
  geom_line()
print(pl)

# Prepare covariate table
covars <- clim_dat %>% 
  select(week_date, week_no, Te, Td, RH_pred) %>% 
  arrange(week_date)

# Create pomp model -------------------------------------------------------
PompMod <- CreateMod(covars_df = covars, lin_bool_val = T)

# Set parameters
pomp::coef(PompMod, names(parms)) <- unname(parms)


# Run deterministic simulation ---------------------------------------------------------
# coef(PompMod, "eps") <- 1
# coef(PompMod, "alpha") <- 1 / (1 * 52)
# coef(PompMod, c("e_Te", "e_RH")) <- -0.2
# coef(PompMod, "R0") <- 1.25
rho_mean_val <- unname(coef(PompMod, "rho_mean"))
rho_k_val <- unname(coef(PompMod, "rho_k"))

# Run trajectory and generate observations 
# NB: this can be done via the rmeasure function in pomp
# but it's actually easier to apply the observation model directly, as follows
sim_det <- trajectory(object = PompMod, format = "data.frame")
sim_det <- sim_det %>% 
  mutate(N_sim = S + I + R, 
         .id = 0,
         CC_obs = rnbinom(n = length(CC), mu = rho_mean_val * CC, size = 1 / rho_k_val))

# Data in long format
sim_det_long <- sim_det %>% 
  pivot_longer(cols = -c(".id", "week"), names_to = "state_var", values_to = "value") %>% 
  mutate(state_var = factor(state_var, levels = c("S", "I", "R", "N_sim", "CC", "CC_obs")))

# Plot variables
pl <- ggplot(data = sim_det_long, 
             mapping = aes(x = week / 52, y = value / parms["N"])) + 
  geom_line() + 
  #scale_y_sqrt() +
  scale_x_continuous(breaks = 0:10) + 
  facet_wrap(~ state_var, scales = "free_y") + 
  labs(x = "Year", y = "Proportion", title = "Deterministic simulation")
print(pl)

# Run stochastic simulation -----------------------------------------------
sim_stoch <- simulate(object = PompMod, nsim = 10, format = "data.frame")
sim_stoch <- sim_stoch %>% 
  mutate(N_sim = S + I + R)

# Data in long format
sim_stoch_long <- sim_stoch %>% 
  pivot_longer(cols = -c(".id", "week"), names_to = "state_var", values_to = "value") %>% 
  mutate(state_var = factor(state_var, levels = c("S", "I", "R", "N_sim", "CC", "CC_obs")))

# Plot variables
pl <- ggplot(data = sim_stoch_long, 
             mapping = aes(x = week / 52, y = value / parms["N"], group = .id)) + 
  geom_line(color = "grey") + 
  #scale_y_sqrt() +
  scale_x_continuous(breaks = 0:10) + 
  facet_wrap(~ state_var, scales = "free_y") + 
  labs(x = "Year", y = "Proportion", title = "Stochastic simulations")
print(pl)

# Plot all ----------------------------------------------------------------
# sim_all <- bind_rows(sim_det_long %>% mutate(type = "det", .id = as.numeric(.id)), 
#                      sim_stoch_long %>%  mutate(type = "stoch", .id = as.numeric(.id)))

pl <- ggplot(data = sim_stoch_long %>% filter(state_var %in% c("S", "I", "R", "CC")), 
             mapping = aes(x = week / 52, y = value / parms["N"], group = .id)) + 
  geom_line(color = "grey") + 
  geom_line(data = sim_det_long %>% filter(state_var %in% c("S", "I", "R", "CC")), color = "red") + 
  scale_x_continuous(breaks = 0:10) + 
  facet_wrap(~ state_var, scales = "free_y") + 
  labs(x = "Year", y = "Proportion", title = "Deterministic (red) vs. stochastic (grey) simulations")
print(pl)


#######################################################################################################
# END
#######################################################################################################
