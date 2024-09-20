####################################################################################################
# Run simulations for vignette on confounding
# Key point: Spatial variability in climate can be a confounder of spatial transmission effects 
# Illustrate by running SIR models with climatic data in different locations with different 
# climates, but no spatial diffusion in transmission. Show that if the variability in climates is 
# not accounted for, it could be wrongly concluded that spatial difussion exists.
# Part 1: Simulations and describing the bias (calculated with ncf, speed with ccf and brms). 

# Out: 
# _saved_/_vignette_spatial_bias/covars_l_{coun_name}.rds
# _saved_/_vignette_spatial_bias/sim_{coun_name}.rds
# _saved/_vignette_spatial_bias/01_brmspeed_{coun_name}.rds
# _saved/_vignette_spatial_bias/vfigD_{coun_name}_simulations.rds
# _saved/_vignette_spatial_bias/vfigB_{coun_name}_sncf.rds
# _saved/_vignette_spatial_bias/vfigC_{coun_name}_ccf.rds
# _saved/_vignette_spatial_bias/vfigA_{coun_name}_map.rds

# _figures/_climate/for each location.pdf
# _figures/_simulations/for each location.pdf
####################################################################################################

# Load packages ------------------------------------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
source("f-Pred_RH.R")
Pred_RH <- Vectorize(FUN = Pred_RH)
source("f-CreateClimData.R")
source("f-PlotClimData.R")
source("f-CreateMod.R")
source("f-SimulateMod.R")
library(pomp)
library(glue)
library(ncf)
library(geosphere)

# Set plot theme
theme_set(theme_classic() + theme(strip.background = element_blank()))

# Set directories 
dirs <- list()
dirs$data <- "_data"
dirs$figures <- "_figures"

# Set model parameters -----------------------------------------------------------------------------
parms <- c("mu" = 1 / 80 / 52, # Birth rate 
           "N" = 5e6, # Total population size
           "R0" = 2.5, # Basic reproduction no
           "sigma_beta" = 0, # SD of environmental noise (set to 0 for a deterministic process model)
           "e_Te" = -0.2, # Effect of Te on transmission 
           "e_RH" = -0.2, # Effect of RH on transmission
           "eps" = 1, # Fraction of infections conferring sterilizing immunity (1: SIR, 0: SIS)
           "alpha" = 1 / (2 * 52), # Rate of waning immunity (to obtain yearly epidemics)
           "rho_mean" = 0.1, # Average reporting probability 
           "rho_k" = 0.04) # Reporting over-dispersion

# Load data ----------------------------------------------------------------------------------------
# Choose between Spain and Colombia
coun_name <- "Colombia"

# Load names of the locations
names_l <- unique(filter(readRDS("_data/spat_data.rds"), country == coun_name)$loc)
# Load spatial data 
spat_dat <- filter(readRDS("_data/spat_data.rds"), country == coun_name)
# Load climatic data per location 
clim_dat_l <- sapply(names_l, simplify = F, function(x) {
  CreateClimData(loc_nm = x, n_years = 10) %>% 
    # Calculate seasonal term of transmission rate
    mutate(beta_seas = exp(parms["e_Te"] * (Te_norm - 1)) * exp(parms["e_RH"] * (RH_pred_norm - 1)),
           Te_norm_lag = lag(x = Te_norm, n = 1L, order_by = week_no), 
           RH_pred_norm_lag = lag(x = RH_pred_norm, n = 1L, order_by = week_no))
})
# Climatic data in long format
clim_dat_long_l <- lapply(clim_dat_l, function(x) {
  pivot_longer(x, cols = Te:RH_pred_norm_miss, names_to = "var", values_to = "value")})

# Check climate and prepare covars -----------------------------------------------------------------

# Check correlation RH and pred RH
unlist(lapply(clim_dat_l, function(x) cor(x$RH_pred_norm, x$RH_norm, method = "spearman"))) %>%
  min()
# Min Spain: LEBL, 0.9376685
# Min Colombia: SKUI, 0.964038
# Min Germany: EDHL, 0.9043391

# Plot climate and save all the locations:
pl_climate_l <- sapply(names_l, fun_PlotClimData, simplify = F)
if(coun_name == "Colombia") pl_climate_l$SKBO else if(coun_name == "Spain") pl_climate_l$LEVS
  
# Plot transmission rate in all the locations 
bind_rows(clim_dat_l, .id = "loc") %>%
  ggplot(aes(x = week_no, y = beta_seas, color = loc)) + 
  geom_line() + 
  scale_color_viridis_d()

# Prepare covariate table 
covars_l <- lapply(clim_dat_l, function(x) {
    select(x, week_date, week_no, Te, Td, RH_pred) %>% 
    arrange(week_date)})

# Save covars df 
saveRDS(covars_l, file = glue("_saved/_vignette_spatial_bias/covars_l_{coun_name}.rds"))

# Run simulations for all locations ----------------------------------------------------------------

# Define parameters to change for the simulations 
chg_parms <- data.frame(eps = 1, 
                        alpha = 1 / (2 * 52),
                        R0 = 2.5) %>% 
  mutate(.id = seq_len(nrow(.)))

# Run simulations
simulations_l <- sapply(names_l, fun_SimulateMod, nsim = 100, parms = parms, 
                        chg_parms = chg_parms, simplify = F)

# Plot simulations 
pl_simulations_l <- sapply(names_l, function(x) {
  pl <- ggplot(data = filter(simulations_l[[x]], sim_id == 1) %>%
                 mutate(id_lab = glue("R0 = {R0}, alpha = 1/52*{1/(52*alpha)}")), 
               aes(x = week_no / 52, y = value / N)) + 
    geom_line() + 
    facet_wrap(~ id_lab + state_var, scales = "free_y", 
               ncol = 4, dir = "v", strip.position = "top") + 
    labs(x = "Year", y = "Proportion") + 
    ggtitle(x)
  ggsave(pl, file = glue("_figures/_simulations/Simulations_{coun_name}_{x}.pdf"), height = 10)
  return(pl)
}, simplify = F)

if(coun_name == "Colombia") pl_simulations_l$SKBO else if(coun_name == "Spain") pl_simulations_l$LEVS

# Save simulations
saveRDS(simulations_l, file = glue("_saved/_vignette_spatial_bias/sim_{coun_name}.rds"))

# Sncf tests =======================================================================================

# Select id of the simulation (id)
set_id <- 1
# Select only 1 simulation 
simulations_1 <- bind_rows(lapply(simulations_l, filter, sim_id == 1)) %>%
  filter(.id == set_id, state_var == "CC_obs") %>%
  mutate(inc = value/N) %>% 
  left_join(spat_dat[c("loc","lon","lat","loc_city_name")]) %>%
  arrange(week_no)

# Prepare the data to wide for sncf
sncf_df <- simulations_1 %>%
  select(loc, lon, lat, week_no, inc) %>%
  pivot_wider(names_from = week_no, values_from = inc)

# Estimate the nonparametric (cross-)correlation function
set.seed(1854)
sncf_sim <- Sncf(x = sncf_df$lon, y = sncf_df$lat, z = sncf_df[,4:523], resamp=1000, latlon = TRUE)

# Spatial spread from CCF lag ======================================================================

# Select reference (distance from more distant)
if(coun_name == "Colombia") reference <- "SKRH" else if(coun_name == "Spain") reference <- "LEAS"

# Prepare the data to wide for ccf
ccf_df <- simulations_1 %>% 
  select(loc, week_no, inc) %>%
  pivot_wider(names_from = loc, values_from = inc)

# Estimate the pair-wise (cross-)correlation function from reference 
ccf_sim_l <- list()
for (i in unique(simulations_1$loc)) {
  ccf_sim_l[["lag"]] <- as.numeric(ccf(ccf_df[[reference]], ccf_df[[i]], lag.max = 35)$lag)
  ccf_sim_l[[i]] <- as.numeric(ccf(ccf_df[[reference]], ccf_df[[i]], lag.max = 35)$acf)
}

# Obtain maximum ccf and distance from reference 
ccf_sim_max <- bind_rows(ccf_sim_l) %>%
  pivot_longer(!lag, names_to = "loc", values_to = "ccf") %>%
  # Filter out negative lags
  filter(lag >= 0) %>%
  # Obtain maximum cross-correlation
  group_by(loc) %>%
  filter(abs(ccf) == max(abs(ccf))) %>%
  ungroup() %>%
  # Obtain distance from reference 
  left_join(spat_dat[c("loc","lon","lat", "lon_km", "lat_km", "loc_city_name")]) %>%
  mutate(lat_km = lat_km + abs(min(lat_km)), 
         lon_km = lon_km + abs(min(lon_km))) 

ccf_sim_max <- ccf_sim_max %>%
  # Calculate distance from reference 
  mutate(dis_km = c(distm(select(ccf_sim_max, lon, lat), 
                          filter(ccf_sim_max, loc == reference)[c("lon","lat")], fun = distGeo)/1000))

# Calculate speed ----------------------------------------------------------------------------------

# Run brms model
brmspeed <- pomp::bake(file = glue("_saved/_vignette_spatial_bias/01_brmspeed_{coun_name}.rds"), {
  brms::brm(lag ~ dis_km, data = ccf_sim_max)
})

# Estimates speed
ccf_speed <- data.frame(brmspeed) %>%
  select(b_dis_km) %>%
  mutate(speed = (1/b_dis_km)*4) %>%
  summarise(mean(speed), sd(speed), quantile(speed, 0.025), quantile(speed, 0.975),
            mean(b_dis_km), sd(b_dis_km), quantile(b_dis_km, 0.025), quantile(b_dis_km, 0.975))
print(ccf_speed)

# Predict from model (for plot)
pred_speed <- predict(brmspeed, 
                      newdata = data.frame(
                        dis_km = seq(min(ccf_sim_max$dis_km), max(ccf_sim_max$dis_km), length.out = 2000))) %>%
  as.data.frame() %>%
  mutate(dis_km = seq(min(ccf_sim_max$dis_km), max(ccf_sim_max$dis_km), length.out = 2000))

#  Main figures ====================================================================================

# Plot the simulations 
pl_main1_inc <- simulations_1 %>%
  group_by(.id, loc) %>%
  mutate(loc_lab = glue("{loc_city_name} ({loc})"),
         inc_minmax = (inc - min(inc))/ (max(inc) - min(inc))) %>%
  ggplot(aes(fill = inc_minmax, x = week_no, y = reorder(loc_lab, lat))) + 
  geom_tile() +
  scale_fill_viridis_c("Scaled incidence", labels = function(x) format(round(x, 1), nsmall = 1)) + 
  scale_x_continuous("Time (weeks)", expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  theme(axis.title.y = element_blank(),
        legend.position = "top",
        strip.background = element_blank())

# Plot sncf
pl_main2_sncf <- data.frame(
  x = c(sncf_sim$real$predicted$x),
  y = c(sncf_sim$real$predicted$y),
  ylow = c(sncf_sim$boot$boot.summary$predicted$y["0.025", ]),
  yhig = c(sncf_sim$boot$boot.summary$predicted$y["0.975", ])) %>%
  mutate(ymean = mean(y),
         yhig = case_when(yhig > 1 ~ 1, .default = yhig),
         ylow = case_when(ylow < -0.1 ~ -0.1, .default = ylow)) %>%
  ggplot(aes(x = x , y = y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yhig), alpha = 0.2, fill = "grey20") + 
  geom_line(color = "grey20", linewidth = 1/2) + 
  geom_hline(aes(yintercept = 0), color = "grey") +
  scale_x_continuous("Distance (km)", expand = c(0,0), limits = c(0,1000)) + 
  scale_y_continuous("Synchrony", limits = c(-0.1,1), expand = c(0,0)) +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank())

# Plot model 
pl_main3_ccf <- ccf_sim_max %>%
  ggplot(aes(y = lag, x = dis_km, color = lag)) +
  geom_ribbon(data = pred_speed,
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5), color = "transparent", alpha = 0.2, fill = "grey20") +
  geom_point(size = 3) +
  geom_smooth(data = pred_speed,
              aes(y = Estimate, color = Estimate), color = "grey20", linewidth = 1/2) +
  scale_color_viridis_c(option = "inferno", direction = -1, end = 0.8, begin = 0.2) +
  scale_y_continuous("Lag (week)") +
  scale_x_continuous("Distance (km)", expand = c(0,0)) +
  facet_grid(~glue(" ")) +
  theme(legend.position = "none",
        strip.background = element_blank())

# Plot map  
pl_main4_map <- ggplot() +
  geom_polygon(data = filter(map_data("world"), region == coun_name),
               aes(x=long, y = lat, group = group), fill = "grey90") +
  geom_point(data = ccf_sim_max,
             aes(x = lon, y = lat, color = lag),
             size = 3) +
  scale_color_viridis_c("Lag \n(week)", option = "inferno", direction = -1, end = 0.8, begin = 0.2) +
  scale_x_continuous("Longitude (°)") +
  scale_y_continuous("Latitude (°N)") +
  facet_grid(~glue("{coun_name}")) +
  coord_map() +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        strip.background = element_blank())

# Save figures 
saveRDS(pl_main1_inc, glue("_saved/_vignette_spatial_bias/vfigD_{coun_name}_simulations.rds"))
saveRDS(pl_main2_sncf, glue("_saved/_vignette_spatial_bias/vfigB_{coun_name}_sncf.rds"))
saveRDS(pl_main3_ccf, glue("_saved/_vignette_spatial_bias/vfigC_{coun_name}_ccf.rds"))
saveRDS(pl_main4_map, glue("_saved/_vignette_spatial_bias/vfigA_{coun_name}_map.rds"))
####################################################################################################
# End
####################################################################################################