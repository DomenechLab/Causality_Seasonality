####################################################################################################
# Run simulations for two-location transmission model with spatial diffusion
# Key point: Spatial variability in climate can be a confounder of spatial transmission effects 
# Illustrate by running SIR models with climatic data in different locations with different 
# climates, but no spatial diffusion in transmission. Show that if the variability in climates is 
# not accounted for, it could be wrongly concluded that spatial difussion exists.
# Part 3: Control for confounding - using two-location transmission model with spatial diffusion

# Out: 
# _saved/_vignette_spatial_bias/_profiles/tj1global_{coun_name}_ref{locref}_{loc}.rds
# _saved/_vignette_spatial_bias/vfigF_{coun_name}_{locref}_profile.rds.rds
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

# Set plot theme
theme_set(theme_classic() + theme(strip.background = element_blank()))

# Set directories 
dirs <- list()
dirs$data <- "_data"
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
# Choose between Spain and Colombia
coun_name <- "Colombia"
# Choose 2 locations
if(coun_name == "Spain") locref <- "LEAS" # Select reference Madrid or Gijon
if(coun_name == "Colombia") locref <- "SKRH" # Select reference Bogota or SKRH

# Load spatial data 
spat_dat <- filter(readRDS("_data/spat_data.rds"), country == coun_name)

# Load covars data
covars_all <- readRDS(file = glue("_saved/_vignette_spatial_bias/covars_l_{coun_name}.rds"))
# Prepare covars data (with reference as location 1)
covars_l <- sapply(str_subset(names(covars_all), locref, negate = TRUE), function(x) {
  covars_all[[locref]] %>%
    rename(Te1 = Te, 
           Td1 = Td, 
           RH_pred1 = RH_pred) %>%
    bind_cols(select(covars_all[[x]], Te, Td, RH_pred)) %>%
    rename(Te2 = Te, 
           Td2 = Td, 
           RH_pred2 = RH_pred)
}, simplify = FALSE)

# Load simulations
simulations_l <- readRDS(file = glue("_saved/_vignette_spatial_bias/sim_{coun_name}.rds"))
# Prepare simulations data 
simulations_trajmatch <- lapply(simulations_l, function(x)
  x %>%
    filter(state_var == "CC_obs") %>%
    rename("sim" = sim_id,
           "CC_obs" = value,
           "week" = week_no) %>% 
    select(sim, week, CC_obs) %>% 
    pivot_wider(names_from = "sim", values_from = "CC_obs") %>% 
    arrange(week) %>% 
    select(-week) %>% 
    as.data.frame())

# Bind two locations (with reference as location 1)
sim_l <- sapply(str_subset(names(covars_all), locref, negate = TRUE), function(x) {
    list(CC_obs1 = simulations_trajmatch[[locref]],
         CC_obs2 = simulations_trajmatch[[x]])
}, simplify = FALSE)

# Plot simulations
bind_rows(simulations_l)  %>%
  filter(state_var == "CC_obs") %>%
  filter(sim_id == 1) %>%
  left_join(spat_dat) %>%
  ggplot(aes(x = week_no, y = value, color = reorder(loc, -lat_km), group = sim_id)) + 
  geom_line(alpha = 1) + 
  facet_grid(reorder(loc, -lat_km)~., scales = "free") + 
  scale_color_viridis_d(direction = -1)

# Generate log-lik profiles for spatial diffusion coupling parameter  ==============================

# Select just the first simulation (of 100)
sim_1 <- lapply(sim_l, function(x) list(CC_obs1 = x$CC_obs1[1], CC_obs2 = x$CC_obs2[1]))

# Global search 
prof_1sim_global <- list()
for(i in names(sim_1)) {
  prof_1sim_global[[i]] <- bake(file = glue(
    "_saved/_vignette_spatial_bias/03_profiles/tj1global_{coun_name}_ref{locref}_{i}.rds"), {
      # FitOptim
      sapply(seq(1e-20, 0.002, length.out = 100), function(asum_coup) {
        FitOptim(sim = sim_1[[i]],
                 all_pars_nm = c("R0", "alpha", "rho_mean", "rho_k", "e_Te", "e_RH"),
                 asum_e_RH = -0.2, 
                 asum_e_Te = -0.2, 
                 asum_coup = asum_coup, 
                 covars = covars_l[[i]])
      }, simplify = F) 
  })
  # Give names
  names(prof_1sim_global[[i]]) <- seq(1e-20, 0.002, length.out = 100)
}

# Summarize global search 
prof_1sim_global_summary <- lapply(prof_1sim_global, function(x) {bind_rows(x, .id = "coup")}) %>%
  bind_rows(.id = "loc") %>%
  group_by(loc) %>%
  mutate(ll_max = max(ll),
         coup_max = coup[which.max(ll == ll_max)]) %>%
  distinct(loc, ll_max, coup_max) 
  
# Check results
lapply(prof_1sim_global, function(x) {bind_rows(x, .id = "coup")}) %>%
  bind_rows(.id = "loc") %>%
  group_by(loc) %>%
  mutate(ll_max = max(ll),
         ll_dif = round(ll - ll_max, digits = 2),
         coup_max = coup[which.max(ll == ll_max)]) %>%
  ggplot(aes(y = ll_dif, x = as.numeric(coup), color = loc, group = loc)) +
    geom_line() + 
    geom_vline(aes(xintercept = as.numeric(coup_max))) + 
    scale_color_viridis_d() + 
    facet_grid(loc~.) 
  
# Local search 
prof_1sim_zoom <- list()
for(i in names(sim_1)) {
  prof_1sim_zoom[[i]] <- bake(file = glue(
    "_saved/_vignette_spatial_bias/03_profiles/tj1local_{coun_name}_ref{locref}_{i}.rds"), {
      # Set limits, lower limit must not be < 0
      upperlim <- as.numeric(filter(prof_1sim_global_summary, loc == i)[["coup_max"]]) + 0.0005
      lowerlim <- as.numeric(filter(prof_1sim_global_summary, loc == i)[["coup_max"]]) - 0.0005
      if(lowerlim < 0) lowerlim <- 1e-20
      # FitOptim
      sapply(seq(lowerlim, upperlim, length.out = 100), function(asum_coup) {
                   FitOptim(sim = sim_1[[i]],
                            all_pars_nm = c("R0", "alpha", "rho_mean", "rho_k", "e_Te", "e_RH"),
                            asum_e_RH = -0.2, 
                            asum_e_Te = -0.2, 
                            asum_coup = asum_coup, 
                            covars = covars_l[[i]])
                 }, simplify = F) 
  })
  # Give names
  upperlim <- as.numeric(filter(prof_1sim_global_summary, loc == i)[["coup_max"]]) + 0.0005
  lowerlim <- as.numeric(filter(prof_1sim_global_summary, loc == i)[["coup_max"]]) - 0.0005
  if(lowerlim < 0) lowerlim <- 1e-20
  names(prof_1sim_zoom[[i]]) <- seq(lowerlim, upperlim, length.out = 100)
}

# Identify CI95%
df_1sim <- bind_rows(bind_rows(lapply(prof_1sim_global, function(x) bind_rows(x, .id = "p")), .id = "loc"), 
                     bind_rows(lapply(prof_1sim_zoom, function(x) bind_rows(x, .id = "p")), .id = "loc")) %>%
  group_by(loc) %>%
  mutate(ll_max = max(ll),
         ll_dif = round(ll - ll_max, digits = 2),
         p_max = p[which.max(ll == ll_max)]) %>%
  ungroup() %>%
  filter(ll_dif >= -2 & ll_dif <= -1.8) %>%
  mutate(p_round = round(as.numeric(p), digits = 4)) %>%
  distinct(loc, p_max, p_round, ll_dif) %>%
  distinct(loc, p_max, p_round) %>%
  group_by(loc) %>%
  mutate(n()) %>%
  filter(p_round == max(p_round))

# Main figure 
pl_main6_profile <- df_1sim %>%
  left_join(spat_dat) %>%
  mutate(loc_lab = glue("{loc_city_name} ({loc})")) %>%
  ggplot(aes(x = as.numeric(p_max), y = reorder(loc_lab, lat_km))) + 
  geom_vline(aes(xintercept = 0), color = "darkorange", linewidth = 1/2) + 
  geom_linerange(aes(xmin = -0.00003, xmax = as.numeric(p_round)), 
                 linewidth = 3, alpha = 0.3, color = "grey20") +
  geom_point(color = "grey20") + 
  scale_y_discrete("") +
  scale_x_continuous(expression(Parameter~tau~estimate), expand = c(0.0001, 0.0001)) + 
  theme(legend.position = "none")

# Save figure
saveRDS(pl_main6_profile, glue("_saved/_vignette_spatial_bias/vfigF_{coun_name}_{locref}_profile.rds"))

####################################################################################################
# END
####################################################################################################
