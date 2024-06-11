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
source("f-CreateClimData.R")
source("f-PlotClimData.R")
source("f-CreateMod.R")
source("f-SimulateMod.R")
library(pomp)
library(glue)
library(ncf)
library(geosphere)

theme_set(theme_classic() + theme(strip.background = element_blank()))

set.seed(1854)

# Set directories 
dirs <- list()
dirs$data <- "_data"
dirs$outputs <- "_outputs"
dirs$figures <- "_figures"

# Set model parameters -----------------------------------------------------------------------------
parms <- c("mu" = 1 / 80 / 52, # Birth rate 
           "N" = 5e6, # Total population size
           "R0" = 2.5, # Basic reproduction no
           "sigma_beta" = 0, # SD of environmental noise (set to 0 for a deterministic process model)
           "e_Te" = -0.2, # Effect of Te on transmission 
           "e_RH" = -0.2, # Effect of RH on transmission
           "eps" = 1, # Fraction of infections conferring sterilizing immunity (1: SIR, 0: SIS)
           "alpha" = 1 / (2 * 52), # Rate of waning immunity
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
pl_climate_l$SKBO
  
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

pl_simulations_l$SKBO
pl_simulations_l$SKPE
# Save simulations
saveRDS(simulations_l, file = glue("_saved/_vignette_spatial_bias/sim_{coun_name}.rds"))

# Simulations Sncf tests ===========================================================================

# Select id of the simulation (id)
set_id = 1
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
sncf_sim <- Sncf(x = sncf_df$lon, y = sncf_df$lat, z = sncf_df[,4:523], resamp=1000, latlon = TRUE)

# Spatial spread from CCF lag ======================================================================

# Select reference (distance from more distant)
if(coun_name == "Colombia") reference <- "SKRH"
if(coun_name == "Spain") reference <- "LEZL"

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

# run brms model
brmspeed <- pomp::bake(file = glue("_saved/_vignette_spatial_bias/brmspeed_{coun_name}.rds"), {
  brms::brm(lag ~ dis_km, data = ccf_sim_max)
})

# estimates speed
ccf_speed <- data.frame(brmspeed) %>%
  select(b_dis_km) %>%
  mutate(speed = (1/b_dis_km)*4) %>%
  summarise(mean(speed), sd(speed), quantile(speed, 0.025), quantile(speed, 0.975),
            mean(b_dis_km), sd(b_dis_km), quantile(b_dis_km, 0.025), quantile(b_dis_km, 0.975))

# predict from model 
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

# Read gp and trajm
pl_gp <- readRDS(glue("_saved/_vignette_spatial_bias/fig_gp.rds"))
pl_gph <- readRDS(glue("_saved/_vignette_spatial_bias/fig_gph_covariance.rds"))
pl_trajm_SKRH <- readRDS(glue("_saved/_vignette_spatial_bias/fig_ci_SKRH.rds")) + 
  ggtitle("Estimation with two-location transmission model") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 1, size = unit(11, "lines")))
pl_trajm_SKBO <- readRDS(glue("_saved/_vignette_spatial_bias/fig_prof_SKBO.rds"))

# # Plots together 
# pl_sim_sncf <- (pl_main4_map| pl_main1_inc + 
#                   theme(legend.justification = "center")) / (pl_main3_ccf | pl_main2_sncf) + 
#   plot_annotation(tag_levels = "A")
# pl_sim_sncf <- (pl_main4_map| pl_main1_inc +
#                   theme(legend.justification = "center")) / (pl_main3_ccf | pl_main2_sncf) / (pl_gp|pl_trajm_SKRH) + 
#   plot_annotation(tag_levels = "A")
                  
pl_up <- ((pl_main4_map| (pl_main2_sncf / pl_main3_ccf) | (pl_main1_inc + theme(legend.justification = "center"))) + 
  plot_layout(widths = c(1,0.5,2)))

pl_vig <- pl_up / (pl_gph + pl_trajm_SKRH + plot_layout(widths = c(2,1), nrow = 1))+ 
  plot_annotation(tag_levels = "A")
print(pl_vig)

# Save figures 
if(coun_name == "Spain") { 
  ggsave(pl_vig, height = 8, width = 10, 
         file = glue("_figures/00_Suppfig_vignette_spatial_bias_{coun_name}.pdf")) 
  }

if(coun_name == "Colombia") { 
  ggsave(pl_vig, height = 8, width = 10, 
         file = glue("_figures/00_Fig_vignette_spatial_bias_{coun_name}.pdf")) 
}

####################################################################################################
# END
####################################################################################################

