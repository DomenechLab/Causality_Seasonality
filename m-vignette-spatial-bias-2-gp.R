####################################################################################################
# Run simulations for vignette on confounding gaussian process models 
# Key point: Climate can confound the estimation of spatial diffusion. 
# We have simulations of populations without spatial diffusion. 
# We want to estimate spatial diffusion using a gaussian process model.
####################################################################################################

# Load packages ------------------------------------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
library(pomp)
library(glue)
library(brms)
library(tidybayes)

theme_set(theme_classic() + theme(strip.background = element_blank()))

set.seed(1854)

# Set directories 
dirs <- list()
dirs$data <- "_data"
dirs$outputs <- "_outputs"
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
brmgp_true <- pomp::bake(file = glue("_saved/_vignette_spatial_bias/brmgp_{coun_name}_true.rds"), {
  brm(CC_obs ~ log(S_lag) + log(I_lag) + Te_norm_lag + RH_pred_norm_lag +
        gp(lat_km, lon_km, scale = FALSE),
      data = df_gp, chains = 4, cores = 4, family = negbinomial(link = "log"),
      iter = 10000, control = list(adapt_delta = 0.9))
})
# Smooth model without climate
brmgp_sconf <- pomp::bake(file = glue("_saved/_vignette_spatial_bias/brmgp_{coun_name}_sconf.rds"), {
  brm(CC_obs ~ s(week_no, k = 50) +
        gp(lat_km, lon_km, scale = FALSE),
      data = df_gp, chains = 4, cores = 4, family = negbinomial(link = "log"),
      iter = 10000, control = list(adapt_delta = 0.9))
})
# Smooth model with climate
brmgp_sclim <- pomp::bake(file = glue("_saved/_vignette_spatial_bias/brmgp_{coun_name}_sclim.rds"), {
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
pl <- post_dfpl %>%
  ggplot(aes(x = x*1000, y = covariance, color = mod_lab)) +
  geom_ribbon(aes(ymin = ci_025_covariance, ymax = ci_975_covariance, fill = mod_lab),
              alpha = 0.1, color = "transparent") + 
  geom_line(aes(group = .draw), linewidth = 1/4, alpha = 0.1) +
  geom_line(aes(y = median_covariance), linewidth = 1/2) +
  scale_x_continuous("Distance (km)", expand = c(0, 0)) +
  scale_y_continuous("Covariance", expand = c(0.0005, 0.0005), limits = c(0,0.04)) +
  scale_color_manual(values = c("darkorange", "grey20", "grey20")) + 
  scale_fill_manual(values = c("darkorange", "grey20", "grey20")) + 
  facet_wrap(.~mod_lab, ncol = 3, strip.position = "top") + 
  ggtitle("Estimation with Gaussian process models") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = unit(11, "lines")))

saveRDS(pl, "_saved/_vignette_spatial_bias/fig_gph_covariance.rds")

# Transform into correlation matrix ----------------------------------------------------------------

# Function to extract correlation matrixes 
fun_gp_cov2cor <- function(etasq = NULL, rhosq = NULL) {
  
  # Number of locations 
  nloc <- length(spat_dist[,1])
  # Compute posterior median covariance among locations
  k <- matrix(0, nrow = nloc, ncol = nloc)
  for (i in 1:nloc) {
    for (j in 1:nloc) {
      k[i, j] <- etasq * exp(-1 * rhosq * spat_dist[i, j]^2)
    }}
  #k <- round(k, 4)
  diag(k) <- etasq + 0.0005
  
  # Convert to correlation matrix
  rho <- round(cov2cor(k), 2)
  colnames(rho) <- rownames(rho) <- colnames(spat_dist)
  
  # Create df with correlations and distances 
  spat_cor <- as.data.frame(rho) %>%
    rownames_to_column("loc") %>%
    pivot_longer(-loc, names_to = "loc2", values_to = "correlation") %>%
    mutate(correlation = replace_na(correlation, 0)) %>%
    left_join(spat_dist %>%
                as.data.frame() %>%
                rownames_to_column("loc") %>%
                pivot_longer(-loc, names_to = "loc2", values_to = "distance"))
  
  return(spat_cor)
}

# Convert to correlations for median and sample of draws -------------------------------------------

# Obtain distance matrix 
spat_dist <- as.matrix(dist(spat_dat[c("lat_km","lon_km")]))/1000
colnames(spat_dist) <- unique(spat_dat$loc)
rownames(spat_dist) <- colnames(spat_dist)

# Correlations for median
med_est <- post_dfpl %>% 
  distinct(mod, median_etasq, median_rhosq, 
           ci_025_etasq, ci_975_etasq, ci_025_rhosq, ci_975_rhosq)

cor_med_l <- list()
cor_med_l <- sapply(c("true", "sconf", "sclim"), function(m_name) {
  fun_gp_cov2cor(etasq = filter(med_est, mod == m_name)[["median_etasq"]],
                 rhosq = filter(med_est, mod == m_name)[["median_rhosq"]])
}, simplify = FALSE) %>%
  bind_rows(.id = "mod") %>%
  mutate(mod_lab = factor(mod, levels = names(mod_labs), labels = mod_labs))

cor_ci_025_l <- list()
cor_ci_025_l <- sapply(c("true", "sconf", "sclim"), function(m_name) {
  fun_gp_cov2cor(etasq = filter(med_est, mod == m_name)[["ci_025_etasq"]],
                 rhosq = filter(med_est, mod == m_name)[["ci_025_rhosq"]])
}, simplify = FALSE) %>%
  bind_rows(.id = "mod") %>%
  mutate(mod_lab = factor(mod, levels = names(mod_labs), labels = mod_labs)) %>%
  rename(correlation_ci_025 = correlation)

cor_ci_975_l <- list()
cor_ci_975_l <- sapply(c("true", "sconf", "sclim"), function(m_name) {
  fun_gp_cov2cor(etasq = filter(med_est, mod == m_name)[["ci_975_etasq"]],
                 rhosq = filter(med_est, mod == m_name)[["ci_975_rhosq"]])
}, simplify = FALSE) %>%
  bind_rows(.id = "mod") %>%
  mutate(mod_lab = factor(mod, levels = names(mod_labs), labels = mod_labs)) %>%
  rename(correlation_ci_975 = correlation)

cor_medci_l <- left_join(cor_med_l, cor_ci_025_l) %>%
  left_join(cor_ci_975_l) %>%
  mutate(mod_lab = factor(mod, levels = names(mod_labs), labels = mod_labs))

# Correlation for a sample of draws 
cor_draws_l <- sapply(c("true", "sconf", "sclim"), function(m_name) {
  lapply(sample(1:length(post_l[[m_name]]$etasq), 100), function(i) {
    fun_gp_cov2cor(etasq = post_l[[m_name]][[i, "etasq"]],
                   rhosq = post_l[[m_name]][[i, "rhosq"]])
  }) %>%
    bind_rows(.id = "draw")
}, simplify = FALSE) %>%
  bind_rows(.id = "mod") %>%
  mutate(mod_lab = factor(mod, levels = names(mod_labs), labels = mod_labs))

# Plot it 
pl_main5_gpcor <- cor_draws_l %>%
  mutate(draw_mod = paste0(mod, draw)) %>%
  ggplot(aes(x = distance * 1000, y = correlation, color = mod_lab, group = draw_mod)) + 
  geom_line(linewidth = 1/4, alpha = 0.1) +
  # geom_ribbon(data = cor_medci_l, 
  #             aes(ymin = correlation_ci_025, ymax = correlation_ci_975,
  #                 fill = mod_lab, group = mod_lab, x = distance * 1000),
  #             alpha = 0.2, color = "transparent") +
  geom_line(data = cor_medci_l, 
            aes(group = mod_lab), linewidth = 1/2) +
  scale_x_continuous("Distance (km)", expand = c(0, 0), limits = c(0, 1000)) +
  scale_y_continuous("Correlation", expand = c(0, 0)) +
  scale_color_brewer(palette = "Set2") + 
  scale_fill_brewer(palette = "Set2") + 
  facet_wrap(.~mod_lab, ncol = 3, strip.position = "top") + 
  theme(legend.position = "none")

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
         `Est.Error` = round(`Est.Error`, digits = 4),
         )

#saveRDS(pl_main5_gpcor, "_saved/_vignette_spatial_bias/fig_gp.rds")
saveRDS(pl_main5_gpcor, "_saved/_vignette_spatial_bias/fig_gph.rds")
####################################################################################################
# END
####################################################################################################

# Correlation matrix
# Number of locations 
etasq <- med_est$median_etasq[2]
rhosq <- med_est$median_rhosq[2]

nloc <- length(spat_dist[,1])
# Compute posterior median covariance among locations
k <- matrix(0, nrow = nloc, ncol = nloc)
for (i in 1:nloc) {
  for (j in 1:nloc) {
    k[i, j] <- etasq * exp(-1 * rhosq * spat_dist[i, j]^2)
  }}
#k <- round(k, 4)
diag(k) <- etasq + 0.0005

# Convert to correlation matrix
rho <- round(cov2cor(k), 2)
colnames(rho) <- rownames(rho) <- colnames(spat_dist)

corrplot::corrplot(rho)

rho["SKBO",]

as.data.frame(spat_dist) %>%
  pivot_longer(everything()) %>%
  filter(value != 0) %>%
  arrange(value)


# Tests 

# Smooth model without climate (DO NOT WORK)
# brmsre_sconf <- pomp::bake(file = glue("_saved/_vignette_spatial_bias/brmsre_{coun_name}_sconf.rds"), {
#   brm(CC_obs ~ s(week_no, k = 50) +
#         s(lat_km, lon_km, bs = "re"),
#       data = df_gp, chains = 4, cores = 4, family = negbinomial(link = "log"),
#       iter = 10000, control = list(adapt_delta = 0.9))
# })
# # Smooth model with climate (DO NOT WORK)
# brmsre_sclim <- pomp::bake(file = glue("_saved/_vignette_spatial_bias/brmsre_{coun_name}_sclim.rds"), {
#   brm(CC_obs ~ s(week_no, k = 50) + Te_norm_lag + RH_pred_norm_lag +
#         s(lat_km, lon_km, bs = "re"),
#       data = df_gp, chains = 4, cores = 4, family = negbinomial(link = "log"),
#       iter = 10000, control = list(adapt_delta = 0.9))
# })

# # Smooth model without climate
# brmre_sconf <- pomp::bake(file = glue("_saved/_vignette_spatial_bias/brmre_{coun_name}_sconf.rds"), {
#   brm(CC_obs ~ s(week_no, k = 50) +
#         (1 + week_no|loc),
#       data = df_gp, chains = 4, cores = 4, family = negbinomial(link = "log"),
#       iter = 10000, control = list(adapt_delta = 0.9))
# })
# 
# # Smooth model with climate
# brmre_sclim <- pomp::bake(file = glue("_saved/_vignette_spatial_bias/brmre_{coun_name}_sclim.rds"), {
#   brm(CC_obs ~ s(week_no, k = 50) + Te_norm_lag + RH_pred_norm_lag +
#         (1 + week_no|loc),
#       data = df_gp, chains = 4, cores = 4, family = negbinomial(link = "log"),
#       iter = 10000, control = list(adapt_delta = 0.9))
# })
# 
# # Smooth model without climate
# brmgpre_sconf <- pomp::bake(file = glue("_saved/_vignette_spatial_bias/brmgpre_{coun_name}_sconf.rds"), {
#   brm(CC_obs ~ s(week_no, k = 50) +
#         gp(lat_km, lon_km, scale = FALSE) + (1 + week_no|loc),
#       data = df_gp, chains = 4, cores = 4, family = negbinomial(link = "log"),
#       iter = 10000, control = list(adapt_delta = 0.9))
# })
# # Smooth model with climate
# brmgpre_sclim <- pomp::bake(file = glue("_saved/_vignette_spatial_bias/brmgpre_{coun_name}_sclim.rds"), {
#   brm(CC_obs ~ s(week_no, k = 50) + Te_norm_lag + RH_pred_norm_lag +
#         gp(lat_km, lon_km, scale = FALSE) + (1 + week_no|loc),
#       data = df_gp, chains = 4, cores = 4, family = negbinomial(link = "log"),
#       iter = 10000, control = list(adapt_delta = 0.9))
# })
