####################################################################################################
# Run simulations for vignette with extra-demographic noise, with 1 or 2 climatic variable forcing
# Key point: Climate can confound the estimation of spatial diffusion. 
####################################################################################################

# Load packages ------------------------------------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
source("f-CreateMod2coup.R")
#source("f-ParticleFilter.R")
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

# Choose between Spain and Colombia
coun_name <- "Colombia"

# Choose 2 locations
if(coun_name == "Spain") covarsnam_2loc <- c("LEVS", "LEAS") # Select only Madrid and Gijon
if(coun_name == "Colombia") covarsnam_2loc <- c("SKBO", "SKPS") # Select only Bogota and Pasto
#if(coun_name == "Colombia") covarsnam_2loc <- c("SKBO", "SKAR") # Select only Bogota and Armenia
#if(coun_name == "Colombia") covarsnam_2loc <- c("SKBO", "SKPE") # Select only Bogota and Pereira

# Load covars data
covars_l <- readRDS(file = glue("_saved/_vignette_spatial_bias/covars_l_{coun_name}.rds"))

# Prepare covariate table
covars <- covars_l[[covarsnam_2loc[1]]] %>%
  rename(Te1 = Te, 
         Td1 = Td, 
         RH_pred1 = RH_pred) %>%
  bind_cols(select(covars_l[[covarsnam_2loc[2]]], Te, Td, RH_pred)) %>%
  rename(Te2 = Te, 
         Td2 = Td, 
         RH_pred2 = RH_pred)

# Simulate with coupling and without extra-demographic noise =======================================

# Create model
PompMod <- CreateMod(covars = covars)
# Define parameters
coef(PompMod)[names(parms)] <- parms

# Define grid of paramaters to change
chg_parms <- expand_grid(p = c(0.001, 0.002, 0.01, 0.02, 0.1, 0.2, 1)) %>% 
  mutate(.id = as.character(seq_len(nrow(.))))

# Create mat of paramaters with changes 
p_mat <- parmat(params = coef(PompMod), nrep = nrow(chg_parms))
p_mat["p", ] <- chg_parms$p

# 100 simulations
sim <- simulate(object = PompMod, params = p_mat, nsim = 100, seed = 2186L, format = "data.frame") %>%
  separate(.id, into = c(".id", "sim")) %>%
  left_join(chg_parms)

# Plot simulations 
sim %>%
  pivot_longer(cols = c("CC_obs1", "CC_obs2")) %>%
  ggplot(aes(x = week, y = value)) + 
  geom_line() +
  facet_grid(name~p) + 
  theme(legend.position = "none") + 
  sim %>%
  filter(sim == "1") %>%
  pivot_longer(cols = c("CC_obs1", "CC_obs2")) %>%
  ggplot(aes(x = week, y = value, color = name)) + 
  geom_line() +
  facet_grid(~p) + 
  theme(legend.position = "none") +
  plot_layout(ncol = 1)

# Save simulations 
saveRDS(filter(sim, sim == 1), file = glue(
  "_saved/_vignette_spatial_bias/sim2c_{coun_name}_{covarsnam_2loc[1]}-{covarsnam_2loc[2]}.rds"))


# Simulate without coupling and with extra-demographic noise =======================================

# Create model
PompMod <- CreateMod(covars = covars)
# Define parameters
coef(PompMod)[names(parms)] <- parms

# Define grid of paramaters to change
if(coun_name == "Colombia") {
  chg_parms <- expand_grid(e_Te = c(0,-0.2)) %>% 
    mutate(.id = as.character(seq_len(nrow(.))))
}
if(coun_name == "Spain") {
  chg_parms <- expand_grid(e_RH = c(0,-0.2)) %>% 
    mutate(.id = as.character(seq_len(nrow(.))))
}

# Create mat of paramaters with changes 
p_mat <- parmat(params = coef(PompMod), nrep = nrow(chg_parms))
if(coun_name == "Colombia") p_mat["e_Te", ] <- chg_parms$e_Te
if(coun_name == "Spain") p_mat["e_RH", ] <- chg_parms$e_RH
  
# 10 simulations
sim <- simulate(object = PompMod, params = p_mat, nsim = 10, seed = 2186L, format = "data.frame") %>%
  separate(.id, into = c(".id", "sim")) %>%
  left_join(chg_parms)

# Plot simulations 
sim %>%
  pivot_longer(cols = c("CC_obs1", "CC_obs2")) %>%
  ggplot(aes(x = week, y = value)) + 
  geom_line() +
  facet_grid(name~.id) + 
  theme(legend.position = "none") + 
sim %>%
  filter(sim == "1") %>%
  pivot_longer(cols = c("CC_obs1", "CC_obs2")) %>%
  ggplot(aes(x = week, y = value, color = name)) + 
  geom_line() +
  facet_grid(~.id) + 
  theme(legend.position = "none") +
plot_layout(ncol = 1)

# Save simulations 
saveRDS(filter(sim, .id == 1), file = glue(
  "_saved/_vignette_spatial_bias/sim2c_{coun_name}_{covarsnam_2loc[1]}-{covarsnam_2loc[2]}_{names(chg_parms)[1]}{chg_parms[1,1]}.rds"))

saveRDS(filter(sim, .id == 2), file = glue(
  "_saved/_vignette_spatial_bias/sim2c_{coun_name}_{covarsnam_2loc[1]}-{covarsnam_2loc[2]}_{names(chg_parms)[1]}{chg_parms[2,1]}.rds"))

####################################################################################################
# END
####################################################################################################



# # Particle filter p and sigma_beta =================================================================
# 
# # P filter assuming -0.2 climate
# test <- fun_ParticleFilter(sim1 = filter(sim, .id == 1), clim = "clim", covars = covars,
#          p_min = 0, p_max = 0.01, p_by = 0.001, 
#          sb_min = 0.02, sb_max = 0.02, sb_by = 1)
# 
# # P filter assuming 0 climate
# test0 <- fun_ParticleFilter(sim1 = filter(sim, .id == 2), clim = "clim0", covars = covars,
#                   p_min = 0.03, p_max = 0.05, p_by = 0.0001, 
#                   sb_min = 0.02, sb_max = 0.02, sb_by = 1)
# 
# # Plot the grid 
# fig_heat_p <- data.frame(ll = c(sapply(test$pfil_l, logLik))) %>%
#   mutate(.id = as.character(1:length(test$pfil_l))) %>%
#   left_join(test$all_parms) %>%
#   mutate(ll_max = max(ll),
#          sb_max = sigma_beta[ll == ll_max],
#          p_max = p[ll == ll_max],
#          clim = "Assuming climate = -0.2") %>%
#   bind_rows(data.frame(ll = c(sapply(test0$pfil_l, logLik))) %>%
#               mutate(.id = as.character(1:length(test0$pfil_l))) %>%
#               left_join(test0$all_parms) %>%
#               mutate(ll_max = max(ll),
#                      sb_max = sigma_beta[ll == ll_max],
#                      p_max = p[ll == ll_max],
#                      clim = "Assuming climate = 0")) %>%
#   filter(p <= 0.1, sigma_beta <= 0.1) %>%
#   ggplot(aes(x = sigma_beta, y = p, fill = ll)) + 
#   geom_tile() + 
#   annotate("rect", xmin=c(0.015), xmax=c(0.025), ymin=c(-0.005), ymax=c(0.005), 
#            colour="white", fill="transparent", size=1) + 
#   geom_point(aes(x = sb_max, y = p_max, group = 1), color = "white") + 
#   scale_y_continuous(expand = c(0,0)) + 
#   scale_x_continuous(expand = c(0,0)) + 
#   scale_fill_viridis_c() + 
#   facet_grid(~clim) + 
#   theme(legend.position = "top")
# 
# fig_prof_p <- data.frame(ll = c(sapply(test$pfil_l, logLik))) %>%
#   mutate(.id = as.character(1:length(test$pfil_l))) %>%
#   left_join(test$all_parms) %>%
#   mutate(ll_max = max(ll),
#          sb_max = sigma_beta[ll == ll_max],
#          p_max = p[ll == ll_max],
#          clim = "Assuming climate = -0.2") %>%
#    bind_rows(data.frame(ll = c(sapply(test0$pfil_l, logLik))) %>%
#               mutate(.id = as.character(1:length(test0$pfil_l))) %>%
#               left_join(test0$all_parms) %>%
#               mutate(ll_max = max(ll),
#                      sb_max = sigma_beta[ll == ll_max],
#                      p_max = p[ll == ll_max],
#                      clim = "Assuming climate = 0")) %>%
#   filter(p <= 0.02) %>%
#   ggplot(aes(x = p, y = ll, color = sigma_beta)) + 
#   geom_line() + 
#   geom_vline(aes(xintercept = 0), color = "grey", linetype = "dotted") + 
#   scale_y_continuous("Profile log-likelihood", expand = c(0,0)) + 
#   scale_x_continuous("Estimate of coupling effect") + 
#   #scale_color_brewer("", palette = "Set2") + 
#   facet_grid(sigma_beta~clim) + 
#   theme(legend.position = "top")
# 
# data.frame(ll = c(sapply(test$pfil_l, logLik))) %>%
#   mutate(.id = as.character(1:length(test$pfil_l))) %>%
#   left_join(test$all_parms) %>%
#   mutate(ll_max = max(ll),
#          sb_max = sigma_beta[ll == ll_max],
#          p_max = p[ll == ll_max],
#          clim = "Assuming climate = -0.2")  %>%
#   arrange(-ll)
# 
# data.frame(ll = c(sapply(test0$pfil_l, logLik))) %>%
#   mutate(.id = as.character(1:length(test0$pfil_l))) %>%
#   left_join(test0$all_parms) %>%
#   mutate(ll_max = max(ll),
#          sb_max = sigma_beta[ll == ll_max],
#          p_max = p[ll == ll_max],
#          clim = "Assuming climate = -0.2")  %>%
#   arrange(-ll)
# 
# data.frame(ll = c(sapply(test0$pfil_l, logLik))) %>%
#   mutate(.id = as.character(1:length(test0$pfil_l))) %>%
#   left_join(test0$all_parms) %>%
#   mutate(ll_max = max(ll),
#          sb_max = sigma_beta[ll == ll_max],
#          p_max = p[ll == ll_max],
#          clim = "Assuming climate = -0.2")  %>%
#   ggplot(aes(x = p, y = ll, group = sigma_beta, color = ll)) + 
#   geom_line() + 
#   geom_vline(aes(xintercept = 0), color = "grey", linetype = "dotted") + 
#   scale_y_continuous("Profile log-likelihood", expand = c(0,0)) + 
#   scale_x_continuous("Estimate of coupling effect") + 
#   facet_wrap(~sigma_beta, scales = "free") + 
#   theme(legend.position = "top")

# test <- `pfil_Colombia_clim_one_p0-1_sb0-1`
# test0 <- `pfil_Colombia_clim0_one_p0-1_sb0-1`
# 
# test <- `pfil_Colombia_clim_one_p0-0.2_sb0-0.2`
# test0 <- `pfil_Colombia_clim0_one_p0-0.2_sb0-0.2`
# 
# test <- `pfil_Colombia_clim_one_p0-1_sb0.06-0.06`
# test0 <- `pfil_Colombia_clim0_one_p0-1_sb0.06-0.06`
# 
# 
# test <- `pfil_Colombia_clim_both_p0-1_sb0-1`
# test0 <- `pfil_Colombia_clim0_both_p0-1_sb0-1`
# 
# test <- `pfil_Colombia_clim_both_p0-0.2_sb0-0.2`
# test0 <- `pfil_Colombia_clim0_both_p0-0.2_sb0-0.2`
# 
# test <- `pfil_Colombia_clim_both_p0-1_sb0.06-0.06`
# test0 <- `pfil_Colombia_clim0_both_p0-1_sb0.06-0.06`
# 
# test <- `pfil_Colombia_clim_p0-0.2_sb0-0.2`
# test0 <- `pfil_Colombia_clim0_p0-0.2_sb0-0.2`
# 
# test <- `pfil_Colombia_clim_p0-1_sb0.02-0.06`
# test0 <- `pfil_Colombia_clim0_p0-1_sb0.02-0.06`

# saveRDS(fig_heat_p, "_saved/_vignette_spatial_bias/fig_heat_p.rds")
# saveRDS(fig_prof_p, "_saved/_vignette_spatial_bias/fig_prof_p.rds")

####################################################################################################
# END
####################################################################################################






































