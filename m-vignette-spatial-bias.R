####################################################################################################
# Run simulations for vignette on confounding
# Key point: spatial variability in climate can be a confounder of spatial
# transmission effects 
# Illustrate by running SIS and SIR models with climatic data in different
# locations with different climates, but no spatial diffusion in transmission.
# Show that if the variability in climates is not accounted for, it could be 
# wrongly concluded that spatial difussion exists (calculated with ncf). 
####################################################################################################

# Load packages ------------------------------------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
source("f-Pred_RH.R")
Pred_RH <- Vectorize(FUN = Pred_RH)
source("f-CreateMod_gau.R")
source("f-CreateClimData.R")
source("f-PlotClimData.R")
source("f-SimulationsMod.R")
source("f-ClimSncf.R")
source("f-SimulationsSncf.R")
source("f-SimulationsSncf_lag.R")

library(pomp)
library(patchwork)
library(glue)
library(ncf)
theme_set(theme_classic())

set.seed(1854)

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
           "rho_mean" = 0.1) # Average reporting probability 
           #"rho_k" = 0.04) # Reporting over-dispersion

# Load names of the locations ----------------------------------------------------------------------
# Choose between Spain and Colombia
country_name <- "Colombia"
names_l <- unique(filter(readRDS("_data/spat_data.rds"), country == country_name)$loc)
#names_l <- str_subset(names_l, "GCLP|Villavicencio|Cartagena|Armenia", negate = T)

# Load climatic data in all the locations ----------------------------------------------------------
e_Te <- parms["e_Te"]
e_RH <- parms["e_RH"]

# Get the dataframes per location
clim_dat_l <- sapply(names_l, simplify = F, function(x) {
  CreateClimData(loc_nm = x, n_years = 10) %>% 
    # Calculate seasonal term of transmission rate
    mutate(beta_seas = exp(e_Te * (Te_norm - 1)) * exp(e_RH * (RH_pred_norm - 1)),
           Te_norm_lag = lag(x = Te_norm, n = 1L, order_by = week_no), 
           RH_pred_norm_lag = lag(x = RH_pred_norm, n = 1L, order_by = week_no))
})

# Check correlation between rh and rh predicted
unlist(lapply(clim_dat_l, function(x) cor(x$RH_pred_norm, x$RH_norm, method = "spearman"))) %>%
  min()
# Min Spain: LEBL, 0.9376685
# Min Colombia: SKUI, 0.964038
# Min Germany: EDHL, 0.9043391

# Data in long format
clim_dat_long_l <- lapply(clim_dat_l, function(x) {
  pivot_longer(x, cols = Te:RH_pred_norm_miss, names_to = "var", values_to = "value")})

# Plots climate ------------------------------------------------------------------------------------
fun_PlotClimData()
# Plot and save all the locations:
# pl_climate_l <- sapply(names_l, fun_PlotClimData, simplify = F)

# Prepare covariate table --------------------------------------------------------------------------
covars_l <- lapply(clim_dat_l, function(x) {
  x %>% 
    select(week_date, week_no, Te, Td, RH_pred) %>% 
    arrange(week_date)})

# Run simulations for all locations ----------------------------------------------------------------

# Define parameters for the simulations 
all_parms <- data.frame(eps = 1, 
                        alpha = 1 / (2 * 52),
                        R0 = 2.5) %>% 
  mutate(.id = seq_len(nrow(.)))

# Run simulations
simulations_l <- sapply(names_l, fun_SimulationsMod, nsim = 100, all_parms = all_parms, simplify = F)

# Plot simulations 
pl_simulations_l <- sapply(names_l, function(x) {
  pl <- ggplot(data = simulations_l[[x]] %>%
                 filter(sim_id == 1) %>%
                 mutate(id_lab = glue("R0 = {R0}, alpha = 1/52*{1/(52*alpha)}")), 
               aes(x = week_no / 52, y = value / N)) + 
    geom_line() + 
    facet_wrap(~ id_lab + state_var, scales = "free_y", 
               ncol = 4, dir = "v", strip.position = "top") + 
    labs(x = "Year", y = "Proportion") + 
    ggtitle(x)
  ggsave(pl, file = glue("_figures/Simulations_{country_name}_{x}.pdf"), height = 10)
  return(pl)
}, simplify = F)

pl_simulations_l$SKBO

# Save simulations
saveRDS(simulations_l, file = glue("_saved/_vignette_spatial_bias/sim_{country_name}.rds"))

# Load spatial data  -------------------------------------------------------------------------------
spat_dat <- readRDS("_data/spat_data.rds") %>%
  filter(country == country_name)

# Climate data Sncf tests --------------------------------------------------------------------------
sncf_clim_l <- list()
sncf_clim_l$Te <- fun_ClimSncf(clim = "Te")
sncf_clim_l$RH_pred <- fun_ClimSncf(clim = "RH_pred")

pl_clim_sncf <- ((sncf_clim_l$Te$plot[[1]] | sncf_clim_l$Te$plot[[2]]) + plot_layout(widths = c(5, 1.5))) / 
  ((sncf_clim_l$RH_pred$plot[[1]] | sncf_clim_l$RH_pred$plot[[2]]) + plot_layout(widths = c(5, 1.5)))
pl_clim_sncf

ggsave(pl_clim_sncf, file = glue("_figures/Sncfclim_{country_name}.pdf"), 
       height = 9, width = 12)

# Simulations Sncf tests ---------------------------------------------------------------------------
sncf_sim <- fun_SimulationsSncf(data_sim_l = lapply(simulations_l, filter, sim_id == 1), set_id = 1)

# Spatial spread from ccf lag ----------------------------------------------------------------------

# Calculate distance from more distant
if(country_name == "Colombia") reference <- "SKRH"
if(country_name == "Spain") reference <- "LEZL"

sncf_sim_lag <- fun_SimulationsSncf_lag(data_sim_l = lapply(simulations_l, filter, sim_id == 1),
                                        set_id = 1, reference = reference)
sncf_sim_lag$speed$`mean(speed)`
sncf_sim_lag$speed$`quantile(speed, 0.025)`
sncf_sim_lag$speed$`quantile(speed, 0.975)`

pl_ss3 <- sncf_sim_lag$data %>%
  ggplot(aes(y = lag, x = dis_km, color = lag)) +
  geom_point(size = 3) +
  geom_smooth(data = sncf_sim_lag$pred,
            aes(y = Estimate, color = Estimate), color = "grey", linewidth = 0.7) +
  geom_ribbon(data = sncf_sim_lag$pred,
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5), color = "transparent", alpha = 0.1) +
  #geom_line(aes(y = lag, x = lag * (sncf_sim_lag$speed$`mean(speed)`/4))) +
  scale_color_viridis_c(option = "inferno", direction = -1, end = 0.8, begin = 0.2) +
  scale_y_continuous("Lag (week)") +
  scale_x_continuous("Distance (km)") +
  facet_grid(~glue(" ")) +
  theme(legend.position = "none",
        strip.background = element_blank())

pl_ss4 <- ggplot() +
  geom_polygon(data = filter(map_data("world"), region == country_name),
               aes(x=long, y = lat, group = group), fill = "grey90") +
  geom_point(data = sncf_sim_lag$data,
             aes(x = lon, y = lat, color = lag),
             size = 3) +
  scale_color_viridis_c("Lag \n(week)", option = "inferno", direction = -1, end = 0.8, begin = 0.2) +
  scale_x_continuous("Longitude (°)") +
  scale_y_continuous("Latitude (°N)") +
  facet_grid(~glue("{country_name}")) +
  coord_map() +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        strip.background = element_blank())

# Plot it! -----------------------------------------------------------------------------------------

pl_sim_sncf <- (pl_ss4| sncf_sim$plot[[1]] + 
                  theme(legend.justification = "center")) / (pl_ss3 | sncf_sim$plot[[2]]) + 
  plot_annotation(tag_levels = "A")


#0.024
if(country_name == "Spain") { 
  ggsave(pl_sim_sncf, height = 7.5, width = 7, 
         file = glue("_figures/00_Suppfig_spatial_diffusion_{country_name}.pdf")) 
  }

if(country_name == "Colombia") { 
  ggsave(pl_sim_sncf, height = 7.5, width = 7, 
         file = glue("_figures/00_Fig_spatial_diffusion_{country_name}.pdf")) 
  }

####################################################################################################
# END
####################################################################################################

