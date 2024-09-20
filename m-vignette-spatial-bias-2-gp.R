####################################################################################################
# Run simulations for vignette on confounding
# Key point: Spatial variability in climate can be a confounder of spatial transmission effects 
# Illustrate by running SIR models with climatic data in different locations with different 
# climates, but no spatial diffusion in transmission. Show that if the variability in climates is 
# not accounted for, it could be wrongly concluded that spatial difussion exists.
# Part 2: Control for confounding - using GP

# Out:
# _saved/_vignette_spatial_bias/02_brmgp_{coun_name}_true.rds
# _saved/_vignette_spatial_bias/02_brmgp_{coun_name}_sconf.rds
# _saved/_vignette_spatial_bias/02_brmgp_{coun_name}_sclim.rds
# _saved/_vignette_spatial_bias/vfigE_{coun_name}_gpcovariance.rds
####################################################################################################

# Load packages ------------------------------------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
library(pomp)
library(glue)
library(brms)
library(tidybayes)

# Set plot theme
theme_set(theme_classic() + theme(strip.background = element_blank()))

# Set directories 
dirs <- list()
dirs$data <- "_data"
dirs$figures <- "_figures"

# Load data ----------------------------------------------------------------------------------------
# Choose between Spain and Colombia
coun_name <- "Colombia"
# Load simulations
simulations_l <- readRDS(file = glue("_saved/_vignette_spatial_bias/sim_{coun_name}.rds"))
# Load spatial data 
spat_dat <- filter(readRDS("_data/spat_data.rds"), country == coun_name)

# Plot them 
bind_rows(lapply(simulations_l, filter, sim_id == 1)) %>%
  filter(state_var == "CC_obs") %>%
  ggplot(aes(x = week_no, y = value, color = loc)) +
  geom_line()

# GP Models ========================================================================================

# Format data for model 
df_gp <- lapply(simulations_l, function(x) {
  # Select only 1 simulation
  filter(x, sim_id == 1) %>%
    pivot_wider(names_from = state_var, values_from = value) %>%
    select(loc, sim_id, week_no, N, S, I, R, CC, CC_obs, Te_norm, RH_pred_norm) %>%
    mutate(S_lag = lag(S, n = 1, order_by = week_no), 
           I_lag = lag(I, n = 1, order_by = week_no),
           Te_norm_lag = lag(x = Te_norm, n = 1L, order_by = week_no), 
           RH_pred_norm_lag = lag(x = RH_pred_norm, n = 1L, order_by = week_no))
}) %>%
  bind_rows() %>%
  left_join(spat_dat)

# Run models ---------------------------------------------------------------------------------------

# True model
brmgp_true <- pomp::bake(file = glue("_saved/_vignette_spatial_bias/02_brmgp_{coun_name}_true.rds"), {
  brm(CC_obs ~ log(S_lag) + log(I_lag) + Te_norm_lag + RH_pred_norm_lag +
        gp(lat_km, lon_km, scale = FALSE),
      data = df_gp, chains = 4, cores = 4, family = negbinomial(link = "log"),
      iter = 10000, control = list(adapt_delta = 0.9))
})
# Smooth model without climate
brmgp_sconf <- pomp::bake(file = glue("_saved/_vignette_spatial_bias/02_brmgp_{coun_name}_sconf.rds"), {
  brm(CC_obs ~ s(week_no, k = 50) +
        gp(lat_km, lon_km, scale = FALSE),
      data = df_gp, chains = 4, cores = 4, family = negbinomial(link = "log"),
      iter = 10000, control = list(adapt_delta = 0.9))
})
# Smooth model with climate
brmgp_sclim <- pomp::bake(file = glue("_saved/_vignette_spatial_bias/02_brmgp_{coun_name}_sclim.rds"), {
  brm(CC_obs ~ s(week_no, k = 50) + Te_norm_lag + RH_pred_norm_lag +
        gp(lat_km, lon_km, scale = FALSE),
      data = df_gp, chains = 4, cores = 4, family = negbinomial(link = "log"),
      iter = 10000, control = list(adapt_delta = 0.9))
})

mod_labs <- c("True Model", "Model with smooth term \n without climate", 
              "Model with smooth term \n including climate")
names(mod_labs) <- c("true", "sconf", "sclim")

# Covariance ---------------------------------------------------------------------------------------

# Obtain draws
post_l <- lapply(list(brmgp_true, brmgp_sconf, brmgp_sclim), function(x) {
  x %>%
    spread_draws(sdgp_gplat_kmlon_km,lscale_gplat_kmlon_km) %>%
    mutate(etasq = sdgp_gplat_kmlon_km^2) %>% 
    mutate(rhosq = 1 / (2 * (lscale_gplat_kmlon_km/1000)^2))
})
post_l <- lapply(list(brmgp_true, brmgp_sconf, brmgp_sclim), function(x) {
  x %>%
    spread_draws(sdgp_gplat_kmlon_km,lscale_gplat_kmlon_km) %>%
    mutate(etasq = sdgp_gplat_kmlon_km^2) %>% 
    mutate(rhosq = 1 / (2 * (lscale_gplat_kmlon_km/1000)^2))
})
names(post_l) <- c("true", "sconf", "sclim")

# Extract covariance 
post_dfpl <- lapply(post_l, function(x) {
  x %>%
    mutate(median_etasq = median(sdgp_gplat_kmlon_km)^2,
           median_rhosq = 1 / (2 * median(lscale_gplat_kmlon_km/1000)^2),
           ci_025_etasq = quantile(sdgp_gplat_kmlon_km, probs = 0.025)^2,
           ci_975_etasq = quantile(sdgp_gplat_kmlon_km, probs = 0.975)^2,
           ci_025_rhosq = 1 / (2 * quantile(lscale_gplat_kmlon_km/1000, probs = 0.025)^2),
           ci_975_rhosq = 1 / (2 * quantile(lscale_gplat_kmlon_km/1000, probs = 0.975)^2)
           ) %>%
    slice_sample(n = 100) %>% 
    select(.draw, etasq, rhosq, median_etasq, median_rhosq, 
           ci_025_etasq, ci_975_etasq, ci_025_rhosq, ci_975_rhosq) %>%
    expand_grid(x = seq(from = 0, to = 0.5, by = .001)) %>%
    mutate(covariance = etasq * exp(-1 * rhosq * x^2),
           median_covariance = median_etasq * exp(-1 * median_rhosq * x^2),
           ci_025_covariance = ci_025_etasq * exp(-1 * ci_025_rhosq * x^2),
           ci_975_covariance = ci_975_etasq * exp(-1 * ci_975_rhosq * x^2),
           )
}) %>%
  bind_rows(.id = "mod") %>%
  mutate(mod_lab = factor(mod, levels = names(mod_labs), labels = mod_labs))

# Plot it 
if(coun_name == "Colombia") {
  uppery <- 0.04
  expany <- 0.0005
} else if(coun_name == "Spain") { 
  uppery <- 0.6
  expany <- 0.005}

# Main figure 
pl_main5_gp <- post_dfpl %>%
  ggplot(aes(x = x*1000, y = covariance, color = mod_lab)) +
  geom_ribbon(aes(ymin = ci_025_covariance, ymax = ci_975_covariance, fill = mod_lab),
              alpha = 0.1, color = "transparent") + 
  geom_line(aes(group = .draw), linewidth = 1/4, alpha = 0.1) +
  geom_line(aes(y = median_covariance), linewidth = 1/2) +
  scale_x_continuous("Distance (km)", expand = c(0, 0)) +
  scale_y_continuous("Covariance", expand = c(expany, expany), limits = c(0, uppery)) +
  scale_color_manual(values = c("darkorange", "grey20", "grey20")) + 
  scale_fill_manual(values = c("darkorange", "grey20", "grey20")) + 
  facet_wrap(.~mod_lab, ncol = 3, strip.position = "top") + 
  ggtitle("Estimation with Gaussian process models") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = unit(11, "lines")))

# Save figure
saveRDS(pl_main5_gp, glue("_saved/_vignette_spatial_bias/vfigE_{coun_name}_gpcovariance.rds"))

# Estimates ----------------------------------------------------------------------------------------
est_l <- lapply(list(brmgp_true, brmgp_sconf, brmgp_sclim), function(x) {
  summary(x)$gp %>%
    rownames_to_column("Paramater")
  }) 
names(est_l) <- c("true", "sconf", "sclim")

est_l %>%
  bind_rows(.id = "mod") %>%
  mutate(`l-95% CI` = round(`l-95% CI`, digits = 4),
         `u-95% CI` = round(`u-95% CI`, digits = 4),
         `Est.Error` = round(`Est.Error`, digits = 4))
####################################################################################################
# End
####################################################################################################
