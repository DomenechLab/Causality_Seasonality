####################################################################################################
# Run simulations for vignette on confounding gaussian process models 
# Key point: Climate can confound the estimation of spatial diffusion. 
# We have simulations of populations without spatial diffusion. 
# We want to estimate spatial diffusion using a gaussian process model.
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
# Choose between Spain and Colombia
coun_name <- "Colombia"
# Choose 2 locations
if(coun_name == "Spain") covarsnam_2loc <- c("LEVS") # Select reference Madrid
if(coun_name == "Colombia") covarsnam_2loc <- c("SKBO", "SKPE") # Select reference Bogota or SKRH

# Load spatial data 
spat_dat <- filter(readRDS("_data/spat_data.rds"), country == coun_name)

# Load covars data
covars_all <- readRDS(file = glue("_saved/_vignette_spatial_bias/covars_l_{coun_name}.rds"))
# Prepare covariate table
covars_l <- covars_all[[covarsnam_2loc[1]]] %>%
    rename(Te1 = Te, 
           Td1 = Td, 
           RH_pred1 = RH_pred) %>%
    bind_cols(select(covars_all[[covarsnam_2loc[2]]], Te, Td, RH_pred)) %>%
    rename(Te2 = Te, 
           Td2 = Td, 
           RH_pred2 = RH_pred)

# Load simulations without extra-demographic noise 
simulations_l <- readRDS(file = glue("_saved/_vignette_spatial_bias/sim_{coun_name}.rds")) 
# Bind two locations 0 p 
sim_0 <- data.frame(
  week = filter(simulations_l[[covarsnam_2loc[1]]], sim_id == 1, state_var == "CC_obs")$week_no,
  CC_obs1 = filter(simulations_l[[covarsnam_2loc[1]]], sim_id == 1, state_var == "CC_obs")$value,
  CC_obs2 = filter(simulations_l[[covarsnam_2loc[2]]], sim_id == 1, state_var == "CC_obs")$value) %>%
  mutate(p = 0)

# Load simulations with p and without extra-demographic noise
sim_all <- readRDS(file = glue(
  "_saved/_vignette_spatial_bias/sim2c_{coun_name}_{covarsnam_2loc[1]}-{covarsnam_2loc[2]}.rds")) %>%
  filter(p != 0) %>%
  # Add old 0 simulations 
  bind_rows(sim_0) %>%
  select(p, week, CC_obs1, CC_obs2) %>%
  arrange(p)

sim_CC_obs1 <- sim_all %>%
    rename("CC_obs" = CC_obs1) %>% 
    select(week, CC_obs, p) %>%
    pivot_wider(names_from = "p", values_from = "CC_obs") %>% 
    arrange(week) %>% 
    select(-week) %>% 
    as.data.frame()

sim_CC_obs2 <- sim_all %>%
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

#saveRDS(sim_l, file = glue("_saved/_vignette_spatial_bias/sim2c_{coun_name}_{covarsnam_2loc[1]}-{covarsnam_2loc[2]}_l.rds"))

# Plot simulations
bind_rows(simulations_l)  %>%
  filter(state_var == "CC_obs") %>%
  filter(sim_id == 1) %>%
  left_join(spat_dat) %>%
  ggplot(aes(x = week_no, y = value, color = reorder(loc, -lat_km), group = sim_id)) + 
  geom_line(alpha = 1) + 
  facet_grid(reorder(loc, -lat_km)~., scales = "free") + 
  scale_color_viridis_d(direction = -1)

# FitOptim =========================================================================================

# Estimate for all p value simulations

# Global profile for 1
prof_1sim_global <- list()
for(i in names(sim_l)) {
  prof_1sim_global[[i]] <- bake(file = glue(
    "_saved/_vignette_spatial_bias/_trajmatch/tjprof_{coun_name}_{covarsnam_2loc[1]}-{covarsnam_2loc[2]}_p{i}.rds"), {
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

prof_1sim_global_summary <- lapply(prof_1sim_global, function(x) {bind_rows(x, .id = "clim")}) %>%
  bind_rows(.id = "p") %>%
  group_by(p) %>%
  mutate(ll_max = max(ll),
         clim_max = clim[which.max(ll == ll_max)]) %>%
  distinct(p, ll_max, clim_max) 

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

# Zoom profile for 1
# 0.035
prof_1sim_zoom <- list()
for(i in names(sim_l)) {
  prof_1sim_zoom[[i]] <- bake(file = glue(
    "_saved/_vignette_spatial_bias/_trajmatch/tjprofzoom_{coun_name}_{covarsnam_2loc[1]}-{covarsnam_2loc[2]}_p{i}.rds"), {
      sapply(seq(as.numeric(filter(prof_1sim_global_summary, p == i)[["clim_max"]]) - 0.040, 
                 as.numeric(filter(prof_1sim_global_summary, p == i)[["clim_max"]]) + 0.040, 
                 length.out = 80), function(asum_clim) {
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
                                    length.out = 80)
}

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

saveRDS(pl, glue("_saved/_vignette_spatial_bias/fig_climconf_{coun_name}_{covarsnam_2loc[1]}-{covarsnam_2loc[2]}.rds"))

####################################################################################################
# END
####################################################################################################

# Estimate for all 100 simulations 
# pars_mle_all <- list()
# for(i in names(sim_l)) {
#   pars_mle_all[[i]] <- FitOptim(sim = sim_l[[i]],
#                            all_pars_nm = c("p","R0", "alpha", "rho_mean", "rho_k", "e_Te", "e_RH"),
#                            asum_clim = -0.2, asum_coup = 1e-20, 
#                            covars = covars_l[[i]])
# }
# 
# pl_est <- bind_rows(pars_mle_all, .id = "loc") %>%
#   mutate(mle = case_when(par %in% c("R0", "alpha", "rho_k") ~ exp(mle),
#                          par %in% c("rho_mean", "p") ~ plogis(mle),
#                          .default = mle),
#          true = case_when(par %in% c("R0", "alpha", "rho_k") ~ exp(true),
#                           par %in% c("rho_mean", "p") ~ plogis(true),
#                           .default = true)) %>%
#   filter(par == "p") %>%
# ggplot(aes(x = mle, color = loc)) + 
#   geom_vline(aes(xintercept = true, color = loc), linetype = "dotted") + 
#   geom_density(alpha = 0.1) + 
#   geom_rug(alpha = 0.5) + 
#   facet_grid(loc~., scales = "free") + 
#   scale_color_brewer("", palette = "Set2") + 
#   scale_x_continuous("Parameter p estimate \n 100 simulations", limits = c(0, 0.001)) + 
#   scale_y_continuous("Density") + 
#   theme(legend.position = "none")











# if(coun_name == "Colombia") {sim_l$sim_sb0.02_oneclim <- readRDS(file = glue(
#   "_saved/_vignette_spatial_bias/sim2c_{coun_name}_{covarsnam_2loc[1]}-{covarsnam_2loc[2]}_e_Te0.rds")) %>%
#   filter(sim == 1)
# sim_l$sim_sb0.02_twoclim <- readRDS(file = glue(
#   "_saved/_vignette_spatial_bias/sim2c_{coun_name}_{covarsnam_2loc[1]}-{covarsnam_2loc[2]}_e_Te-0.2.rds")) %>%
#   filter(sim == 1)}
# if(coun_name == "Spain") {sim_l$sim_sb0.02_oneclim <- readRDS(file = glue(
#   "_saved/_vignette_spatial_bias/sim2c_{coun_name}_{covarsnam_2loc[1]}-{covarsnam_2loc[2]}_e_RH0.rds")) %>%
#   filter(sim == 1)
# sim_l$sim_sb0.02_twoclim <- readRDS(file = glue(
#   "_saved/_vignette_spatial_bias/sim2c_{coun_name}_{covarsnam_2loc[1]}-{covarsnam_2loc[2]}_e_RH-0.2.rds")) %>%
#   filter(sim == 1)}












# 
# 
# # Select just the first simulation 
# sim_1 <- lapply(sim_l, function(x) {
#   list(CC_obs1 = x$CC_obs1[1],
#        CC_obs2 = x$CC_obs2[1])
# })
# 
# # IF USING ALL LOCATIONS WITH SKRH 
# if (locref == "SKRH") {
#   # Global profile for 1
#   prof_1sim_global <- list()
#   for(i in names(sim_1)) {
#     prof_1sim_global[[i]] <- bake(file = glue(
#       "_saved/_vignette_spatial_bias/_trajmatch/tj1prof_{coun_name}_ref{locref}_{i}.rds"), {
#         sapply(seq(1e-20, 0.002, length.out = 100), function(asum_coup) {
#           FitOptim(sim = sim_1[[i]],
#                    all_pars_nm = c("R0", "alpha", "rho_mean", "rho_k", "e_Te", "e_RH"),
#                    asum_clim = -0.2, asum_coup = asum_coup, 
#                    covars = covars_l[[i]])
#         }, simplify = F) 
#     })
#     names(prof_1sim_global[[i]]) <- seq(1e-20, 0.002, length.out = 100)
#   }
# 
#   # Zoom profile for 1
#   prof_1sim_zoom <- list()
#   for(i in names(sim_1)) {
#     prof_1sim_zoom[[i]] <- bake(file = glue(
#       "_saved/_vignette_spatial_bias/_trajmatch/tj1profz_{coun_name}_ref{locref}_{i}.rds"), {
#         sapply(seq(0.00020, 0.00159, length.out = 100), function(asum_coup) {
#           FitOptim(sim = sim_1[[i]],
#                    all_pars_nm = c("R0", "alpha", "rho_mean", "rho_k", "e_Te", "e_RH"),
#                    asum_clim = -0.2, asum_coup = asum_coup, 
#                    covars = covars_l[[i]])
#         }, simplify = F) 
#       })
#     names(prof_1sim_zoom[[i]]) <- seq(0.00020, 0.00159, length.out = 100)
#   }
#   
#   # Bind to df
#   df_global <- bind_rows(lapply(prof_1sim_global, function(x) bind_rows(x, .id = "p")), .id = "loc")
#   df_zoom <- bind_rows(lapply(prof_1sim_zoom, function(x) bind_rows(x, .id = "p")), .id = "loc")
#   
#   # Identify CI95%
#   df_1sim <- bind_rows(df_global, df_zoom) %>%
#     group_by(loc) %>%
#     mutate(ll_max = max(ll),
#            ll_dif = round(ll - ll_max, digits = 2),
#            p_max = p[which.max(ll == ll_max)]) %>%
#     ungroup() %>%
#     filter(ll_dif >= -2.0 & ll_dif < -1.84) %>%
#     mutate(p_round = round(as.numeric(p), digits = 4)) %>%
#     filter(!(loc == "SKBO" & ll_dif %in% c(-1.86, -2))) %>%
#     filter(!(loc == "SKCG" & ll_dif %in% c(-1.86, -2))) %>%
#     filter(!(loc == "SKUI" & ll_dif == -1.97)) %>%
#     filter(!(loc == "SKVV" & p_round %in% c(0.0012, 0.0011))) %>%
#     filter(!(loc == "SKMD" & p_round == 0.0006)) %>%
#     filter(!(loc == "SKPE" & p_round == 0.0004)) %>%
#     distinct(loc, p_max, p_round) %>%
#     group_by(loc) %>%
#     mutate(n())
# 
#   pl <- df_1sim %>%
#     left_join(spat_dat) %>%
#     mutate(loc_lab = glue("{loc_city_name} ({loc})")) %>%
#     ggplot(aes(x = as.numeric(p_max), y = reorder(loc_lab, lat_km))) + 
#     geom_vline(aes(xintercept = 0), color = "grey", linetype = "dotted") + 
#     geom_linerange(aes(xmin = 0, xmax = as.numeric(p_round)), 
#                    size = 3, alpha = 0.3, color = "dodgerblue3") +
#     geom_point(color = "dodgerblue3") + 
#     scale_y_discrete("") +
#     scale_x_continuous("Parameter p estimate") + 
#     theme(legend.position = "none")
#   
#   saveRDS(pl, "_saved/_vignette_spatial_bias/fig_ci_SKRH.rds")
# }
# 
# # IF USING SOME LOCATIONS FOR SKBO 
# if (locref == "SKBO") {
#   # Profile for 1
#   prof_1sim <- list()
#   for(i in names(sim_1)) {
#     prof_1sim[[i]] <- bake(file = glue(
#       "_saved/_vignette_spatial_bias/_trajmatch/tj1prof_{coun_name}_ref{locref}_{i}.rds"), {
#         sapply(seq(1e-20, 0.002, length.out = 40), function(asum_coup) {
#           FitOptim(sim = sim_1[[i]],
#                    all_pars_nm = c("R0", "alpha", "rho_mean", "rho_k", "e_Te", "e_RH"),
#                    asum_clim = -0.2, asum_coup = asum_coup, 
#                    covars = covars_l[[i]])
#         }, simplify = F) 
#       })
#     names(prof_1sim[[i]]) <- seq(1e-20, 0.002, length.out = 40)
#   }
#   
#   df_1sim <- bind_rows(lapply(prof_1sim, function(x) bind_rows(x, .id = "p")), .id = "loc") %>%
#     group_by(loc) %>%
#     mutate(ll_max = max(ll),
#            ll_dif = round(ll - ll_max, digits = 2),
#            p_max = p[which.max(ll == ll_max)]) %>%
#     ungroup()
# 
#   pl <- df_1sim %>%
#     left_join(spat_dat) %>%
#     filter(loc %in% c("SKPE", "SKAR", "SKPS")) %>%
#     mutate(loc_lab = glue("{loc_city_name} ({loc})")) %>%
#     ggplot(aes(x = round(as.numeric(p), digits = 5), y = ll_dif, color = reorder(loc_lab, -lat_km))) + 
#     geom_hline(aes(yintercept = -1.92), color = "grey") + 
#     geom_smooth(linewidth = 1/2, se = FALSE) + 
#     scale_color_viridis_d("", end = 0.8) + 
#     scale_y_continuous("Log-likelihood - max log-likehood", limits = c(-4, 0)) + 
#     scale_x_continuous("Parameter p estimate") + 
#     geom_vline(aes(xintercept = as.numeric(p_max), color = reorder(loc_lab, lat_km))) + 
#     theme(legend.position = c(0.9, 0.9))
#   
#   saveRDS(pl, "_saved/_vignette_spatial_bias/fig_prof_SKBO.rds")
# }












