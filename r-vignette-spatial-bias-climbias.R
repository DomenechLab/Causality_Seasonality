####################################################################################################
# Run simulations for two-location transmission model with spatial diffusion
# Key point: Reviewer point: Can spatial diffusion bias the effect of climate on transmission? 
####################################################################################################

# Load packages ------------------------------------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
source("f-CreateMod2coup.R")
source("f-RunTrajMatch.R")
library("bbmle")
library(pomp)
library(glue)
library(patchwork)

theme_set(theme_classic() + theme(strip.background = element_blank()))

# Set directories 
dirs <- list()
dirs$data <- "_data"
dirs$outputs <- "_outputs"
dirs$figures <- "_figures"

# Load climatic data in a given location------------------------------------------------------------
parms <- c("mu" = 1 / 80 / 52, # Birth rate 
           "p" = 0, # Coupling 
           "N1" = 5e6, # Total population size
           "N2" = 5e6, # Total population size
           "R0" = 2.5, # Reproduction no
           "sigma_beta" = 0, # SD of environmental noise
           "e_Te" = -0.2, # Effect of temperature
           "e_RH" = -0.2, # Effect of RH
           "eps" = 1, # Fraction of infections conferring immunity
           "alpha" = 1 / (2 * 52), # Waning rate 
           "rho_mean" = 0.1, # Average reporting probability
           "rho_k" = 0.04) # Reporting over-dispersion

# Load climatic data in a given location -----------------------------------------------------------
# Choose Colombia
coun_name <- "Colombia"
# Choose 2 locations
if(coun_name == "Colombia") covarsnam_2loc <- c("SKBO", "SKPE") # Select reference Bogota or SKRH

# Load spatial data 
spat_dat <- filter(readRDS("_data/spat_data.rds"), country == coun_name)

# Load covars data
covars_all <- readRDS(file = glue("_saved/_vignette_spatial_bias/covars_l_{coun_name}.rds"))
# Prepare covars data
covars_l <- covars_all[[covarsnam_2loc[1]]] %>%
    rename(Te1 = Te, 
           Td1 = Td, 
           RH_pred1 = RH_pred) %>%
    bind_cols(select(covars_all[[covarsnam_2loc[2]]], Te, Td, RH_pred)) %>%
    rename(Te2 = Te, 
           Td2 = Td, 
           RH_pred2 = RH_pred)

# Load simulations without spatial diffusion -------------------------------------------------------
simulations_l <- readRDS(file = glue("_saved/_vignette_spatial_bias/sim_{coun_name}.rds")) 
# Bind two locations 0 p 
sim_0 <- data.frame(
  week = filter(simulations_l[[covarsnam_2loc[1]]], sim_id == 1, state_var == "CC_obs")$week_no,
  CC_obs1 = filter(simulations_l[[covarsnam_2loc[1]]], sim_id == 1, state_var == "CC_obs")$value,
  CC_obs2 = filter(simulations_l[[covarsnam_2loc[2]]], sim_id == 1, state_var == "CC_obs")$value) %>%
  mutate(p = 0)

# Make simulations with spatial diffusion ----------------------------------------------------------

# Create model
PompMod <- CreateMod(covars = covars_l)
# Define parameters
coef(PompMod)[names(parms)] <- parms
# Define grid of paramaters to change
chg_parms <- expand_grid(p = c(0.001, 0.002, 0.01, 0.02, 0.1, 0.2, 1)) %>% 
  mutate(.id = as.character(seq_len(nrow(.))))
# Create mat of paramaters with changes 
p_mat <- parmat(params = coef(PompMod), nrep = nrow(chg_parms))
p_mat["p", ] <- chg_parms$p

# 1 simulation
sim <- simulate(object = PompMod, params = p_mat, nsim = 1, seed = 2186L, format = "data.frame") %>%
  left_join(chg_parms) %>%
  # Add old 0 simulations 
  bind_rows(sim_0) %>%
  select(p, week, CC_obs1, CC_obs2) %>%
  arrange(p)

# Plot simulations 
sim %>%
  pivot_longer(cols = c("CC_obs1", "CC_obs2")) %>%
  ggplot(aes(x = week, y = value, color = name)) + 
  geom_line() +
  facet_grid(p~.) + 
  theme(legend.position = "top") +
  plot_layout(ncol = 1)

# Prepare simulations for back-fitting -------------------------------------------------------------     
sim_CC_obs1 <- sim %>%
    rename("CC_obs" = CC_obs1) %>% 
    select(week, CC_obs, p) %>%
    pivot_wider(names_from = "p", values_from = "CC_obs") %>% 
    arrange(week) %>% 
    select(-week) %>% 
    as.data.frame()

sim_CC_obs2 <- sim %>%
  rename("CC_obs" = CC_obs2) %>% 
  select(week, CC_obs, p) %>%
  pivot_wider(names_from = "p", values_from = "CC_obs") %>% 
  arrange(week) %>% 
  select(-week) %>% 
  as.data.frame()

sim_l <- sapply(colnames(sim_CC_obs1), function(x) {
  list(CC_obs1 = data.frame(sim_CC_obs1[,x]),
       CC_obs2 = data.frame(sim_CC_obs2[,x]))
}, simplify = FALSE)

# Estimate the effect of climate assuming no spatial diffusion for all p value simulations =========

# Global search
prof_1sim_global <- list()
for(i in names(sim_l)) {
  prof_1sim_global[[i]] <- bake(file = glue(
    "_saved/_vignette_spatial_bias/03_profiles/r_climbias/tj1global_{coun_name}_{covarsnam_2loc[1]}-{covarsnam_2loc[2]}_p{i}.rds"), {
      sapply(seq(-0.25, 0, length.out = 40), function(asum_clim) {
        FitOptim(sim = sim_l[[i]],
                 all_pars_nm = c("R0", "alpha", "rho_mean", "rho_k", "e_Te"),
                 asum_e_RH = asum_clim, 
                 asum_e_Te = -0.2,
                 asum_coup = 0, 
                 covars = covars_l)
      }, simplify = F) 
    })
  names(prof_1sim_global[[i]]) <- seq(-0.25, 0, length.out = 40)
}

# Summary of global search
prof_1sim_global_summary <- lapply(prof_1sim_global, function(x) {bind_rows(x, .id = "clim")}) %>%
  bind_rows(.id = "p") %>%
  group_by(p) %>%
  mutate(ll_max = max(ll),
         clim_max = clim[which.max(ll == ll_max)]) %>%
  distinct(p, ll_max, clim_max) 

# Plot check
lapply(prof_1sim_global, function(x) {bind_rows(x, .id = "clim")}) %>%
  bind_rows(.id = "p") %>%
  group_by(p) %>%
  mutate(ll_max = max(ll),
         ll_dif = round(ll - ll_max, digits = 2),
         clim_max = clim[which.max(ll == ll_max)]) %>%
  ggplot(aes(y = ll_dif, x = as.numeric(clim), color = as.numeric(p), group = p)) +
    geom_line() + 
    geom_vline(aes(xintercept = as.numeric(clim_max))) + 
    scale_color_viridis_c() + 
    facet_grid(p~.) 

# Local search
prof_1sim_zoom <- list()
for(i in names(sim_l)) {
  prof_1sim_zoom[[i]] <- bake(file = glue(
    "_saved/_vignette_spatial_bias/03_profiles/r_climbias/tj1local_{coun_name}_{covarsnam_2loc[1]}-{covarsnam_2loc[2]}_p{i}.rds"), {
      sapply(seq(as.numeric(filter(prof_1sim_global_summary, p == i)[["clim_max"]]) - 0.040, 
                 as.numeric(filter(prof_1sim_global_summary, p == i)[["clim_max"]]) + 0.040, 
                 length.out = 100), function(asum_clim) {
        FitOptim(sim = sim_l[[i]],
                 all_pars_nm = c("R0", "alpha", "rho_mean", "rho_k", "e_Te"),
                 asum_e_RH = asum_clim, 
                 asum_e_Te = -0.2,
                 asum_coup = 0, 
                 covars = covars_l)
      }, simplify = F) 
    })
  names(prof_1sim_zoom[[i]]) <- seq(as.numeric(filter(prof_1sim_global_summary, p == i)[["clim_max"]]) - 0.040, 
                                    as.numeric(filter(prof_1sim_global_summary, p == i)[["clim_max"]]) + 0.040, 
                                    length.out = 100)
}

# Figure for the review 
pl <- lapply(prof_1sim_zoom, function(x) {bind_rows(x, .id = "clim")}) %>%
  bind_rows(.id = "p") %>%
  group_by(p) %>%
  mutate(ll_max = max(ll),
         ll_dif = round(ll - ll_max, digits = 2),
         clim_max = clim[which.max(ll == ll_max)]) %>%
  ggplot(aes(y = ll_dif, x = as.numeric(clim), color = p, group = p)) +
  geom_hline(aes(yintercept = -1.92), color = "grey") + 
  geom_vline(aes(xintercept = -0.2), color = "grey") + 
  geom_line() + 
  geom_vline(aes(xintercept = as.numeric(clim_max), color = p)) + 
  scale_color_viridis_d() + 
  scale_x_continuous("Estimated effect of relative humidity") + 
  scale_y_continuous("Profile log-likelihood") + 
  facet_grid(p~.) +
  theme(legend.position = "none") +
  ggtitle(glue("{coun_name}: Bogotá (SKBO) - Pereira (SKPE)"))

print(pl)

####################################################################################################
# END
####################################################################################################
