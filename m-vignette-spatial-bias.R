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
source("f-SimulationsSncf_lag.R")

library(pomp)
library(patchwork)
library(glue)
library(ncf)
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

# Load names of the locations ----------------------------------------------------------------------
country_name <- "Colombia"
names_l <- unique(filter(readRDS("_data/clim_data.rds"), country == country_name)$loc)
names_l <- str_subset(names_l, "GCLP|Villavicencio|Cartagena|Armenia", negate = T)

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
clim_dat_long_l <- lapply(clim_dat_l, function(x) {
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
all_parms <- data.frame(eps = c(1, 1, 1, 1), 
                        alpha = c(1 / (1 * 52), 1 / (1 * 52),  1 / (2 * 52),  1 / (2 * 52)),
                        R0 = c(1.25, 2.5, 1.25, 2.5)) %>% 
  mutate(.id = seq_len(nrow(.)))

# Run simulations
simulations_l <- sapply(names_l, fun_SimulationsMod, all_parms = all_parms, simplify = F)

# Plot simulations 
pl_simulations_l <- sapply(names_l, function(x) {
  pl <- ggplot(data = simulations_l[[x]] %>%
                 mutate(id_lab = glue("SIR, R0 = {R0}, alpha = 1/52*{1/(52*alpha)}")), aes(x = week_no / 52, y = value / N)) + 
    geom_line() + 
    facet_wrap(~ id_lab + state_var, scales = "free_y", 
               ncol = 4, dir = "v", strip.position = "top") + 
    labs(x = "Year", y = "Proportion") + 
    ggtitle(x)
  ggsave(pl, file = glue("_figures/{country_name}_{x}_simulations.pdf"), height = 10)
  return(pl)
}, simplify = F)

pl_simulations_l$SKAR

# Load spatial data  -------------------------------------------------------------------------------
spat_dat <- readRDS("_data/spat_data.rds") %>%
  filter(loc_country_name == country_name)

# Plot climate data organized by latitude 
bind_rows(clim_dat_long_l) %>% 
  left_join(spat_dat[c("loc","lon","lat","loc_city_name")]) %>%
  filter(var %in% c("Te_norm", "Td_norm", "RH_norm", "RH_pred_norm")) %>%
  ggplot(aes(x = week_date, y = value, color = reorder(loc,-lat))) + 
  geom_line() + 
  scale_color_viridis_d("Station", option = "turbo") + 
  facet_grid(reorder(loc,-lat) ~ var, scales = "free")

# Climate data Sncf tests --------------------------------------------------------------------------

sncf_clim_l <- list()
sncf_clim_l$Te <- fun_ClimSncf(clim = "Te")
sncf_clim_l$RH_pred <- fun_ClimSncf(clim = "RH_pred")

pl_clim_sncf <- ((sncf_clim_l$Te$plot[[1]] | sncf_clim_l$Te$plot[[2]]) + plot_layout(widths = c(5, 1.5))) / 
  ((sncf_clim_l$RH_pred$plot[[1]] | sncf_clim_l$RH_pred$plot[[2]]) + plot_layout(widths = c(5, 1.5)))
pl_clim_sncf

ggsave(pl_clim_sncf, file = glue("_figures/{country_name}_sncfclim.pdf"), 
       height = 9, width = 12)

# Simulations Sncf tests ---------------------------------------------------------------------------

sncf_sim_l <- list()
sncf_sim_l$id1 <- fun_SimulationsSncf(sim_id = 1)
sncf_sim_l$id2 <- fun_SimulationsSncf(sim_id = 2)
sncf_sim_l$id3 <- fun_SimulationsSncf(sim_id = 3)
sncf_sim_l$id4 <- fun_SimulationsSncf(sim_id = 4)

sncf_sim_l$id4$plot

# Spatial spread -----------------------------------------------------------------------------------

# Spatial spread from peak timing ------------------------------------------------------------------

# Create spatial spread data.frame
ss_df <- bind_rows(simulations_l) %>%
  filter(state_var == "CC") %>%
  # Select trajectory
  filter(.id == 4) %>%
  # Add spatial data 
  left_join(spat_dat[c("loc","lon","lat","loc_city_name")]) %>%
  mutate(loc_lab = glue("{loc_city_name} ({loc})")) %>%
  group_by(.id, loc) %>%
  # Incidence, week and year
  mutate(inc = value/N,
         week_yr = rep(1:52, 10),
         year = rep(1:10, each = 52)) %>%
  group_by(.id, loc, loc_lab, lat, week_yr) %>%
  mutate(mean_inc = mean(inc)) %>%
  ungroup() %>%
  # Peak for the mean seasonality
  group_by(.id, loc, loc_lab, lat) %>%
  mutate(max_mean_inc = max(mean_inc)) %>%
  ungroup() %>%
  # Peak for each year 
  group_by(.id, loc, loc_lab, lat, year) %>%
  mutate(max_inc = max(inc)) 

# Plot maximum values of incidence 
ss_df %>%
  ggplot(aes(y = inc, x = week_yr, group = year, color = year)) + 
  geom_line() +
  geom_line(aes(x = week_yr, y = max_inc, group = year, color = year)) + 
  scale_color_viridis_c("Year") + 
  scale_x_continuous("Week", expand = c(0,0)) + 
  scale_y_continuous("Incidence", expand = c(0,0)) + 
  facet_wrap(reorder(loc_lab, -lat)~., scales = "free_y") + 
  ggtitle("R0 = 2.5, alpha = 1/2*52")

# Transform degrees to km 
sslm_df <- ss_df %>%
  filter(max_inc == inc) %>%
  distinct(lon, lat, loc_lab, loc, week_yr, year) %>%
  ungroup() %>%
  mutate(lat_km = (dsm::latlong2km(lat = lat, lon = lon))$km.n,
         lat_km = lat_km + abs(min(lat_km)), 
         lon_km = (dsm::latlong2km(lat = lat, lon = lon))$km.e,
         lon_km = lon_km + abs(min(lon_km)))
  
# Calculate distance from more distant 
if(country_name == "Colombia") reference <- "SKRH"
if(country_name == "Spain") reference <- "LEZL"

lim_lat <- mean(unlist(filter(sslm_df, loc == reference)[,"lat_km"]))
lim_lon <- mean(unlist(filter(sslm_df, loc == reference)[,"lon_km"]))
if(country_name == "Spain") sslm_df <- filter(sslm_df, year != 1)

sslm_df <- sslm_df %>%
  mutate(dis_km = sqrt((lat_km - lim_lat)^2+(lon_km  - lim_lon)^2)) 

mod <- lm(week_yr~dis_km, data = sslm_df)
speed <- list(peak_speed = list(mean = (1/coef(mod)[[2]])*4,
                                sd = (1/coef(mod)[[2]])*4 - (1/(coef(mod)[[2]] + sqrt(diag(vcov(mod)))[[2]]))*4)) #km/mo

pl_ss1 <- sslm_df %>%
  ggplot(aes(y = dis_km, x = week_yr, color = week_yr)) + 
  geom_point(size = 2) + 
  scale_color_viridis_c(option = "inferno", direction = -1, end = 0.8, begin = 0.2) + 
  scale_x_continuous("Peak timing (week)") + 
  scale_y_continuous("Distance (km)")+ 
  facet_grid(~glue("Speed mean = {round(speed$peak_speed$mean)}, std. error = {round(speed$peak_speed$sd)} km/mo")) + 
  theme(legend.position = "none",
        strip.background = element_blank())

pl_ss2 <- ggplot() +
  geom_polygon(data = filter(map_data("world"), region == country_name), 
               aes(x=long, y = lat, group = group), fill = "grey90") + 
  geom_point(data = sslm_df %>%
               group_by(lon, lat, loc_lab) %>%
               summarise(mean_week_yr = mean(week_yr)),
             aes(x = lon, y = lat, color = mean_week_yr),
             size = 2) + 
  scale_color_viridis_c("Peak timing \n(week)", option = "inferno", direction = -1, end = 0.8, begin = 0.2) + 
  scale_x_continuous("Longitude (°)") + 
  scale_y_continuous("Latitude (°N)") + 
  facet_grid(~glue("{country_name}")) +
  coord_map() + 
  theme(panel.grid = element_blank(),
        legend.position = "left",
        strip.background = element_blank())

# Spatial spread from ccf lag ----------------------------------------------------------------------

sncf_sim_lag <- fun_SimulationsSncf_lag(sim_id = 4, reference = reference)
speed$ccf_speed <- sncf_sim_lag$speed

pl_ss3 <- sncf_sim_lag$data %>%
  ggplot(aes(y = dis_km, x = lag, color = lag)) + 
  geom_point(size = 2) + 
  scale_color_viridis_c(option = "inferno", direction = -1, end = 0.8, begin = 0.2) + 
  scale_x_continuous("Lag (week)") + 
  scale_y_continuous("Distance (km)") + 
  facet_grid(~glue("Speed mean = {round(speed$ccf_speed$mean)}, std. error = {round(speed$ccf_speed$sd)} km/mo")) + 
  theme(legend.position = "none",
        strip.background = element_blank())

pl_ss4 <- ggplot() +
  geom_polygon(data = filter(map_data("world"), region == country_name), 
               aes(x=long, y = lat, group = group), fill = "grey90") + 
  geom_point(data = sncf_sim_lag$data,
             aes(x = lon, y = lat, color = lag),
             size = 2) + 
  scale_color_viridis_c("Lag \n(week)", option = "inferno", direction = -1, end = 0.8, begin = 0.2) + 
  scale_x_continuous("Longitude (°)") + 
  scale_y_continuous("Latitude (°N)") + 
  facet_grid(~glue("{country_name}")) +
  coord_map() + 
  theme(panel.grid = element_blank(),
        legend.position = "left",
        strip.background = element_blank())

# Plot it! -----------------------------------------------------------------------------------------

pl_sim_sncf <- ((sncf_sim_l$id4$plot[[1]] | sncf_sim_l$id4$plot[[2]]) + plot_layout(widths = c(3,1))) / (pl_ss1 | pl_ss2 | pl_ss3 | pl_ss4) 
pl_sim_sncf <- sncf_sim_l$id4$plot[[1]] / (sncf_sim_l$id4$plot[[2]] | pl_ss3 | pl_ss4 )

ggsave(pl_sim_sncf, file = glue("_figures/{country_name}_sncfsimulations.pdf"), 
       height = 9, width = 15)

####################################################################################################
# END
####################################################################################################



