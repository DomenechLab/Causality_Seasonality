####################################################################################################
# Run simulations for vignette on the total causal effect of temperature 
# Key point: As temperature affects humidity, the total causal effect of temperature comprises 
# direct and indirect effects (mediated by humidity).
# Illustration: 
####################################################################################################

# Load packages ------------------------------------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
source("f-Pred_RH.R")
Pred_RH <- Vectorize(FUN = Pred_RH)
source("f-CreateMod.R")
source("f-CreateClimData.R")
source("f-PlotClimData.R")
source("f-SimulationsMod.R")

library(pomp)
library(patchwork)
library(glue)
theme_set(theme_classic())

# Set directories 
dirs <- list()
dirs$data <- "_data"
dirs$outputs <- "_outputs"
dirs$figures <- "_figures"

# Set model parameters -----------------------------------------------------------------------------
parms <- c("mu" = 1 / 80 / 52, # Birth rate 
           "N" = 5e6, # Total population size
           "R0" = 1.25, # Basic reproduction no
           "e_Te" = -0.2, # Effect of Te on transmission 
           "e_RH" = -0.2, # Effect of RH on transmission
           "eps" = 1, # Fraction of infections conferring sterilizing immunity (1: SIR, 0: SIS)
           "alpha" = 0, # Rate of waning immunity
           "rho_mean" = 0.1, # Average reporting probability 
           "rho_k" = 0.04) # Reporting over-dispersion

# Load climatic data in a given location------------------------------------------------------------
# Bogota (Colombia): weather station "SKBO"
# Rostock (Germany): Rostock
loc_nm <- "Rostock" ## Choose here the location!

if(loc_nm == "SKBO") loc_tidy_nm <- "Bogotá"
if(loc_nm == "Rostock") loc_tidy_nm <- "Rostock"

clim_dat <- CreateClimData(loc_nm = loc_nm, n_years = 10) 

e_Te <- parms["e_Te"]
e_RH <- parms["e_RH"]

# Calculate seasonal term of transmission rate
clim_dat <- clim_dat %>% 
  mutate(beta_seas = exp(e_Te * (Te_norm - 1)) * exp(e_RH * (RH_pred_norm - 1)))

# Data in long format
clim_dat_long <- clim_dat %>% 
  pivot_longer(cols = Te:RH_pred_norm_miss, names_to = "var", values_to = "value")

# Plot time series of climatic variables
fun_PlotClimData(loc_nm, clim_data_list = F, country = "Germany")

# Add lagged variables for Te_norm and RH_norm
# clim_dat <- clim_dat %>% 
#   mutate(Te_norm_lag = lag(x = Te_norm, n = 1L, order_by = week_no), 
#          RH_pred_norm_lag = lag(x = RH_pred_norm, n = 1L, order_by = week_no))

# Prepare covariate table --------------------------------------------------------------------------
covars <- clim_dat %>% 
  select(week_date, week_no, Te, Td, RH_pred) %>% 
  arrange(week_date)

# Create pomp model --------------------------------------------------------------------------------
PompMod <- CreateMod(covars = covars, lin_bool_val = T)

# Set parameters
pomp::coef(PompMod, names(parms)) <- unname(parms)
base_pars <- pomp::coef(PompMod)

# Run simulations ----------------------------------------------------------------------------------
pomp::coef(PompMod, names(base_pars)) <- unname(base_pars)
rho_mean_val <- unname(coef(PompMod, "rho_mean"))
rho_k_val <- 0.04 # Reporting overdispersion
N_val <- unname(coef(PompMod, "N"))

R0_seq <- c(1.25, 2.5)
alpha_seq <- 1 / 52
eps_seq <- 1
e_Te_seq <- c(0, -0.2)
e_RH_seq <- -0.2

all_parms <- expand_grid(R0 = R0_seq, alpha = alpha_seq, eps = eps_seq, e_Te = e_Te_seq, e_RH = e_RH_seq) %>% 
  mutate(.id = seq_len(nrow(.)))

p_mat <- parmat(params = coef(PompMod), 
                nrep = nrow(all_parms))
p_mat["eps", ] <- all_parms$eps
p_mat["R0", ] <- all_parms$R0
p_mat["alpha", ] <- all_parms$alpha
p_mat["e_Te", ] <- all_parms$e_Te
p_mat["e_RH", ] <- all_parms$e_RH

# Run trajectory and generate observations ---------------------------------------------------------
sim_raw <- trajectory(object = PompMod,
                      params = p_mat,
                      format = "data.frame")

sim <- sim_raw %>% 
  mutate(.id = as.integer(.id)) %>% 
  left_join(y = all_parms) %>% 
  mutate(N_sim = S + I + R, 
         CC_obs = rnbinom(n = length(CC), mu = rho_mean_val * CC, size = 1 / rho_k_val)) %>% 
  rename(week_no = week) %>% 
  left_join(y = clim_dat %>% select(week_no, matches("Te_norm|RH_pred_norm"), beta_seas, Te, RH)) %>% 
  select(.id, eps, R0, alpha, week_no, S, I, R, N_sim, CC, CC_obs, everything())

# Data in long format
sim_long <- sim %>% 
  pivot_longer(cols = c(S:CC_obs), names_to = "state_var", values_to = "value") %>% 
  mutate(state_var = factor(state_var, levels = c("S", "I", "R", "N_sim", "CC", "CC_obs", "beta_seas")), 
         N = N_val) %>%
  mutate(Te_effect = case_when((e_Te == 0 & e_RH == -0.2) ~ "Indirect effect",
                               (e_Te == -0.2 & e_RH == 0) ~ "Direct effect",
                               (e_Te == -0.2 & e_RH == -0.2) ~ "Total effect",
                               (e_Te == 0 & e_RH == 0) ~ "No effect")) %>%
  mutate(te_rh = factor(glue("Te = {e_Te}, RH = {e_RH}"))) %>%
  mutate(Te_effect = glue("{Te_effect} ({te_rh})")) %>%
  filter(te_rh != "Te = 0, RH = 0") %>% 
  mutate(beta_seas = exp(e_Te * (Te_norm - 1)) * exp(e_RH * (RH_pred_norm - 1))) 

# Plot simulations
if(loc_nm == "SKBO") col <- c("#8CD1BB", "#407967")
if(loc_nm == "Rostock") col <- c("pink", "firebrick2")

sim_long %>%
  ggplot(aes(x = week_no / 52, y = value / N, color = Te_effect, group = Te_effect)) + 
  geom_line() + 
  facet_grid(state_var~R0, scales = "free_y") + 
  labs(x = "Year", y = "Proportion") + 
  scale_color_manual("Effect of temperature", values = col) + 
  theme(legend.position = "top")

# Plot effect of temperature -----------------------------------------------------------------------
pl1 <- sim_long %>%
  filter(state_var == "CC") %>%
  filter(R0 == 1.25) %>%
  ggplot(aes(x = week_no, y = value / N * 100, color = Te_effect, group = Te_effect)) + 
  geom_line() + 
  scale_x_continuous("Time (weeks)") + 
  scale_y_continuous("Total cases (per week per 100)") + 
  scale_color_manual("Effect of temperature", values = col) + 
  theme(strip.background = element_blank(),
        legend.position = "none")

pl2 <- sim_long %>%
  filter(R0 == 1.25) %>% 
  ggplot(aes(x = week_no, y = beta_seas, color = Te_effect, group = Te_effect)) + 
  geom_line() + 
  scale_x_continuous("Time (weeks)") + 
  scale_y_continuous("Transmission rate") + 
  scale_color_manual("Effect of temperature", values = col) +
  facet_grid(~glue("{loc_tidy_nm}")) + 
  theme(legend.position = "top",
        legend.direction = "vertical",
        strip.background = element_blank())

pl3 <- sim_long %>%
  filter(state_var == "CC") %>%
  filter(R0 == 1.25) %>% 
  #  mutate(Te_cut = cut(Te, breaks = 20, labels = FALSE),
  #        Te_cut_lab = cut(Te, breaks = 20)) %>%
  # group_by(.id, Te_effect, Te_cut, Te_cut_lab) %>%
  # summarise(Te = mean(Te),
  #           CC_min = min(value/N*100),
  #           CC_max = max(value/N*100),
  #           CC = mean(value/N*100)) %>%
  ggplot(aes(x = Te, y = value/N*100, color = Te_effect, group = Te_effect)) + 
  geom_point() + 
  #geom_linerange(aes(ymin = CC_min, ymax = CC_max)) + 
  scale_x_continuous("Temperature") + 
  scale_y_continuous("Total cases (per week per 100)") + 
  facet_wrap(~Te_effect, scales = "free") + 
  scale_color_manual("Effect of temperature", values = col) + 
  geom_smooth(color = "grey") + 
  theme(legend.position = "none",
        strip.background = element_blank())

pl4 <- sim_long %>%
  filter(R0 == 1.25) %>% 
  mutate(Te_cut = cut(Te, breaks = 10, labels = FALSE),
         Te_cut_lab = cut(Te, breaks = 10)) %>%
  ggplot(aes(x = Te, y = beta_seas, color = Te_effect, group = Te_effect)) + 
  geom_point() + 
  scale_x_continuous("Temperature") + 
  scale_y_continuous("Transmission rate") + 
  facet_wrap(~Te_effect, scales = "free") + 
  scale_color_manual("Effect of temperature", values = col) + 
  geom_smooth(color = "grey") + 
  theme(legend.position = "none",
        strip.background = element_blank())

pl_tot <- (pl2 + pl4) / (pl1 + pl3) 
pl_tot

ggsave(pl_tot, file = glue("_figures/vignette-total-effect-{loc_tidy_nm}.pdf"), 
       height = 7, width = 12)



