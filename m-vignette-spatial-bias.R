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
source("f-CreateMod.R")
source("f-CreateClimData.R")
source("f-PlotClimData.R")
source("f-SimulationsMod.R")
source("f-ClimSncf.R")
source("f-SimulationsSncf.R")

library(pomp)
library(patchwork)
library(glue)
library(ncf)
theme_set(theme_bw())

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

# Load names of the locations ----------------------------------------------------------------------
country_name <- "Colombia"
names_l <- unique(filter(readRDS("_data/clim_data.rds"), country == country_name)$loc)
names_l <- str_subset(names_l, "GCLP|Villavicencio", negate = T)

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

# Data in long format
clim_dat_log_l <- lapply(clim_dat_l, function(x) {
  pivot_longer(x, cols = Te:RH_pred_norm_miss, names_to = "var", values_to = "value")})

# Plots climate ------------------------------------------------------------------------------------
fun_PlotClimData()
pl_climate_l <- sapply(names_l, fun_PlotClimData, simplify = F)


# Prepare covariate table --------------------------------------------------------------------------
covars_l <- lapply(clim_dat_l, function(x) {
  x %>% 
    select(week_date, week_no, Te, Td, RH_pred) %>% 
    arrange(week_date)})

# Run simulations for all locations ----------------------------------------------------------------

# Define parameters for the simulations 
all_parms <- data.frame(eps = c(0, 1, 1, 1, 1), 
                        alpha = c(0, 0, 0,  1 / (1 * 52),  2 / (1 * 52)),
                        R0 = c(1.25, 1.25, 2.5, 2.5, 2.5)) %>% 
  mutate(.id = seq_len(nrow(.)))

# Run simulations
simulations_l <- sapply(names_l, fun_SimulationsMod, all_parms = all_parms, simplify = F)

# Plot simulations 
pl_simulations_l <- sapply(names_l, function(x) {
  ggplot(data = simulations_l[[x]], aes(x = week_no / 52, y = value / N)) + 
    geom_line() + 
    facet_wrap(~ .id + state_var, scales = "free_y", 
               ncol = 5, dir = "v", strip.position = "right") + 
    labs(x = "Year", y = "Proportion") + 
    ggtitle(x)
}, simplify = F)

pl_simulations_l$Armenia

# Load spatial data  -------------------------------------------------------------------------------
spat_dat <- readRDS("_data/spat_data.rds") %>%
  filter(loc_country_name == country_name)

# Plot climate data organized by latitude 
bind_rows(clim_dat_log_l) %>% 
  left_join(spat_dat[c("loc","lon","lat","loc_city_name")]) %>%
  filter(var %in% c("Te_norm", "Td_norm", "RH_norm", "RH_pred_norm")) %>%
  ggplot(aes(x = week_date, y = value, color = reorder(loc,-lon))) + 
  geom_line() + 
  scale_color_viridis_d("Station", option = "turbo") + 
  facet_grid(reorder(loc,-lon) ~ var, scales = "free")

# Climate data Sncf tests --------------------------------------------------------------------------

sncf_clim_l <- list()
sncf_clim_l$Te <- fun_ClimSncf(clim = "Te")
sncf_clim_l$Td <- fun_ClimSncf(clim = "Td")
sncf_clim_l$RH_pred <- fun_ClimSncf(clim = "RH_pred")

pl_clim_sncf <- sncf_clim_l$Te$plot / sncf_clim_l$Td$plot / sncf_clim_l$RH_pred$plot
pl_clim_sncf

ggsave(pl_clim_sncf, file = glue("_figures/{country_name}_sncfclim.pdf"), 
       height = 12, width = 12*3/5)

# Simulations Sncf tests ---------------------------------------------------------------------------

sncf_sim_l <- list()
sncf_sim_l$id1 <- fun_SimulationsSncf(sim_id = 1)
sncf_sim_l$id2 <- fun_SimulationsSncf(sim_id = 2)
sncf_sim_l$id3 <- fun_SimulationsSncf(sim_id = 3)
sncf_sim_l$id4 <- fun_SimulationsSncf(sim_id = 4)
sncf_sim_l$id5 <- fun_SimulationsSncf(sim_id = 5)

pl_sim_sncf <- sncf_sim_l$id1$plot / sncf_sim_l$id2$plot / sncf_sim_l$id3$plot / 
  sncf_sim_l$id4$plot / sncf_sim_l$id5$plot

pl_sim_sncf

ggsave(pl_sim_sncf, file = glue("_figures/{country_name}_sncfsimulations.pdf"), 
       height = 12, width = 12)

####################################################################################################
# END
####################################################################################################























