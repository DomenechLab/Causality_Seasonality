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

# Choose between Spain and Colombia
country_name <- "Colombia"

# Read simulations ---------------------------------------------------------------------------------
simulations_l <- readRDS(file = glue("_saved/_vignette_spatial_bias/sim_{country_name}.rds"))

# Extract the simulations from Bogota
loc1 = simulations_l$SKBO %>%
  pivot_wider(names_from = state_var, values_from = value) %>%
  select(.id, sim_id, week_no, N, S, I, R, CC, CC_obs, Te_norm, RH_pred_norm)
names(loc1) <- paste0(names(loc1), "_loc1")

# Extract the simulations from Pasto
loc2 = simulations_l$SKPS %>%
  pivot_wider(names_from = state_var, values_from = value) %>%
  select(.id, sim_id, week_no, N, S, I, R, CC, CC_obs, Te_norm, RH_pred_norm)
names(loc2) <- paste0(names(loc2), "_loc2")

# Get the dataframes together 
dat_2loc = bind_cols(loc1, loc2) %>%
  filter(sim_id_loc1 == sim_id_loc2,
         week_no_loc1 == week_no_loc2) %>%
  rename(sim = sim_id_loc1,
         week_no = week_no_loc1) %>%
  select(-c(sim_id_loc2, week_no_loc2))

# Plot them 
dat_2loc %>%
  filter(sim == 1) %>%
  ggplot(aes(x = week_no, y = CC_obs_loc1)) +
  geom_line() +
  geom_line(aes(y = CC_obs_loc2), color = "red")

# Run models ---------------------------------------------------------------------------------------
library(mgcv)

M_est <- list()
M_est_se <- list()
R2 <- list()
M_cor <- list()
for(i in 1:100) {
  # Prepare covariates 
  df_mod <- dat_2loc %>%
    filter(sim == i) %>%
    mutate(S_lag_loc1 = lag(S_loc1, n = 1, order_by = week_no), 
           I_lag_loc1 = lag(I_loc1, n = 1, order_by = week_no),
           Te_norm_lag_loc1 = lag(x = Te_norm_loc1, n = 1L, order_by = week_no), 
           RH_pred_norm_lag_loc1 = lag(x = RH_pred_norm_loc1, n = 1L, order_by = week_no),
           S_lag_loc2 = lag(S_loc2, n = 1, order_by = week_no), 
           I_lag_loc2 = lag(I_loc2, n = 1, order_by = week_no),
           Te_norm_lag_loc2 = lag(x = Te_norm_loc2, n = 1L, order_by = week_no), 
           RH_pred_norm_lag_loc2 = lag(x = RH_pred_norm_loc2, n = 1L, order_by = week_no))
  
  # Models 
  mod_true <- gam(list(log(CC_obs_loc1) ~ log(S_lag_loc1) + log(I_lag_loc1) + Te_norm_lag_loc1 + RH_pred_norm_lag_loc1,
                       log(CC_obs_loc2) ~ log(S_lag_loc2) + log(I_lag_loc2) + Te_norm_lag_loc2 + RH_pred_norm_lag_loc2),
             family = mvn(d = 2), data = df_mod)
  mod_conf <- gam(list(log(CC_obs_loc1) ~ log(S_lag_loc1) + log(I_lag_loc1),
                       log(CC_obs_loc2) ~ log(S_lag_loc2) + log(I_lag_loc2)),
             family = mvn(d = 2), data = df_mod)
  mod_scli <- gam(list(log(CC_obs_loc1) ~ s(week_no, k = 50) + Te_norm_lag_loc1 + RH_pred_norm_lag_loc1,
                       log(CC_obs_loc2) ~ s(week_no, k = 50) + Te_norm_lag_loc2 + RH_pred_norm_lag_loc2),
                  family = mvn(d = 2), data = df_mod)
  mod_smoo <- gam(list(log(CC_obs_loc1) ~ s(week_no, k = 50),
                       log(CC_obs_loc2) ~ s(week_no, k = 50)),
                  family = mvn(d = 2), data = df_mod)
  
  # Save estimates 
  M_est[["true"]][[i]] <- coef(mod_true)
  M_est_se[["true"]][[i]] <- summary(mod_true)$se
  R2[["true"]][[i]] <- summary(mod_true)$r.sq
  M_cor[["true"]][[i]] <- cov2cor(solve(crossprod(mod_true$family$data$R)))[1,2]
  
  M_est[["conf"]][[i]] <- coef(mod_conf)
  M_est_se[["conf"]][[i]] <- summary(mod_conf)$se
  R2[["conf"]][[i]] <- summary(mod_conf)$r.sq
  M_cor[["conf"]][[i]] <- cov2cor(solve(crossprod(mod_conf$family$data$R)))[1,2]
 
  M_est[["scli"]][[i]] <- coef(mod_scli)
  M_est_se[["scli"]][[i]] <- summary(mod_scli)$se
  R2[["scli"]][[i]] <- summary(mod_scli)$r.sq
  M_cor[["scli"]][[i]] <- cov2cor(solve(crossprod(mod_scli$family$data$R)))[1,2]
  
  M_est[["smoo"]][[i]] <- coef(mod_smoo)
  M_est_se[["smoo"]][[i]] <- summary(mod_smoo)$se
  R2[["smoo"]][[i]] <- summary(mod_smoo)$r.sq
  M_cor[["smoo"]][[i]] <- cov2cor(solve(crossprod(mod_smoo$family$data$R)))[1,2]
}

as.data.frame(lapply(M_cor, function(x) mean(unlist(x))))
# true      conf      scli      smoo
# 1 -0.01991325 0.6866993 0.6040748 0.6181526

# Obtain confidence intervals ----------------------------------------------------------------------
library(MASS)

# Consider just one simulation 
df_mod <- dat_2loc %>%
  filter(sim == 1) %>%
  mutate(S_lag_loc1 = lag(S_loc1, n = 1, order_by = week_no), 
         I_lag_loc1 = lag(I_loc1, n = 1, order_by = week_no),
         Te_norm_lag_loc1 = lag(x = Te_norm_loc1, n = 1L, order_by = week_no), 
         RH_pred_norm_lag_loc1 = lag(x = RH_pred_norm_loc1, n = 1L, order_by = week_no),
         S_lag_loc2 = lag(S_loc2, n = 1, order_by = week_no), 
         I_lag_loc2 = lag(I_loc2, n = 1, order_by = week_no),
         Te_norm_lag_loc2 = lag(x = Te_norm_loc2, n = 1L, order_by = week_no), 
         RH_pred_norm_lag_loc2 = lag(x = RH_pred_norm_loc2, n = 1L, order_by = week_no))


# Obtain confidence intervals for the TRUE model --------------------------------------------------- 
mod_true <- gam(list(log(CC_obs_loc1) ~ log(S_lag_loc1) + log(I_lag_loc1) + Te_norm_lag_loc1 + RH_pred_norm_lag_loc1,
                     log(CC_obs_loc2) ~ log(S_lag_loc2) + log(I_lag_loc2) + Te_norm_lag_loc2 + RH_pred_norm_lag_loc2),
                family = mvn(d = 2), data = df_mod)
  
# Get Sigma
cov_mat_y = solve(crossprod(mod_true$family$data$R))
# Get the means 
mu_pred = predict(mod_true)
# Simulate 1
ysim = mvrnorm(n = 100, mu = mu_pred[1, ], Sigma = cov_mat_y)
# Simulate all 
ysim <- list()
for (i in 1:length(mu_pred[ ,1])) {
  ysim[[i]] = as.data.frame(mvrnorm(n = 100, mu = mu_pred[i, ], Sigma = cov_mat_y)) %>%
    rownames_to_column(var = "sim")
}
# Bind all simulations 
ysim_df <- bind_rows(ysim, .id = "week_no")
  
corre_l <- list()
for (sim_i in 1:length(unique(ysim_df$sim))) {
  # Add to df 
  df_mod_sim <- df_mod %>%
    mutate(loc1_sim = c(NA, filter(ysim_df, sim == sim_i)$V1),
           loc2_sim = c(NA, filter(ysim_df, sim == sim_i)$V2),
           loc1_sim = exp(loc1_sim),
           loc2_sim = exp(loc2_sim))
  
  # Run the model again 
  gam_model_sim <- gam(list(log(loc1_sim) ~ log(S_lag_loc1) + log(I_lag_loc1) + Te_norm_lag_loc1 + RH_pred_norm_lag_loc1,
                            log(loc2_sim) ~ log(S_lag_loc2) + log(I_lag_loc2) + Te_norm_lag_loc2 + RH_pred_norm_lag_loc2),
                       family = mvn(d = 2), data = df_mod_sim)
  # Extract correlation 
  corre_l[[sim_i]] = cov2cor(solve(crossprod(gam_model_sim$family$data$R)))[1,2]
}

mean(unlist(corre_l))
quantile(unlist(corre_l), c(0.025, 0.975))

# Obtain confidence intervals for the CONFOUNDED model ---------------------------------------------
mod_conf <- gam(list(log(CC_obs_loc1) ~ log(S_lag_loc1) + log(I_lag_loc1),
                     log(CC_obs_loc2) ~ log(S_lag_loc2) + log(I_lag_loc2)),
                family = mvn(d = 2), data = df_mod)

# Get Sigma
cov_mat_y = solve(crossprod(mod_conf$family$data$R))
# Get the means 
mu_pred = predict(mod_conf)
# Simulate 1
ysim = mvrnorm(n = 100, mu = mu_pred[1, ], Sigma = cov_mat_y)
# Simulate all 
ysim <- list()
for (i in 1:length(mu_pred[ ,1])) {
  ysim[[i]] = as.data.frame(mvrnorm(n = 100, mu = mu_pred[i, ], Sigma = cov_mat_y)) %>%
    rownames_to_column(var = "sim")
}
# Bind all simulations 
ysim_df <- bind_rows(ysim, .id = "week_no")

corre_l <- list()
for (sim_i in 1:length(unique(ysim_df$sim))) {
  # Add to df 
  df_mod_sim <- df_mod %>%
    mutate(loc1_sim = c(NA, filter(ysim_df, sim == sim_i)$V1),
           loc2_sim = c(NA, filter(ysim_df, sim == sim_i)$V2),
           loc1_sim = exp(loc1_sim),
           loc2_sim = exp(loc2_sim))
  
  # Run the model again 
  gam_model_sim <- gam(list(log(loc1_sim) ~ log(S_lag_loc1) + log(I_lag_loc1),
                            log(loc2_sim) ~ log(S_lag_loc2) + log(I_lag_loc2)),
                       family = mvn(d = 2), data = df_mod_sim)
  # Extract correlation 
  corre_l[[sim_i]] = cov2cor(solve(crossprod(gam_model_sim$family$data$R)))[1,2]
}

mean(unlist(corre_l))
quantile(unlist(corre_l), c(0.025, 0.975))

# Obtain confidence intervals for the SMOOTH WITH CLIMATE model ------------------------------------

mod_scli <- gam(list(log(CC_obs_loc1) ~ s(week_no, k = 50) + Te_norm_lag_loc1 + RH_pred_norm_lag_loc1,
                     log(CC_obs_loc2) ~ s(week_no, k = 50) + Te_norm_lag_loc2 + RH_pred_norm_lag_loc2),
                family = mvn(d = 2), data = df_mod)

# Get Sigma
cov_mat_y = solve(crossprod(mod_scli$family$data$R))
# Get the means 
mu_pred = predict(mod_scli)
# Simulate 1
ysim = mvrnorm(n = 100, mu = mu_pred[1, ], Sigma = cov_mat_y)
# Simulate all 
ysim <- list()
for (i in 1:length(mu_pred[ ,1])) {
  ysim[[i]] = as.data.frame(mvrnorm(n = 100, mu = mu_pred[i, ], Sigma = cov_mat_y)) %>%
    rownames_to_column(var = "sim")
}
# Bind all simulations 
ysim_df <- bind_rows(ysim, .id = "week_no")

corre_l <- list()
for (sim_i in 1:length(unique(ysim_df$sim))) {
  # Add to df 
  df_mod_sim <- df_mod %>%
    mutate(loc1_sim = c(NA, filter(ysim_df, sim == sim_i)$V1),
           loc2_sim = c(NA, filter(ysim_df, sim == sim_i)$V2),
           loc1_sim = exp(loc1_sim),
           loc2_sim = exp(loc2_sim))
  
  # Run the model again 
  gam_model_sim <- gam(list(log(loc1_sim) ~ s(week_no, k = 50) + Te_norm_lag_loc1 + RH_pred_norm_lag_loc1,
                            log(loc2_sim) ~ s(week_no, k = 50) + Te_norm_lag_loc2 + RH_pred_norm_lag_loc2),
                       family = mvn(d = 2), data = df_mod_sim)
  # Extract correlation 
  corre_l[[sim_i]] = cov2cor(solve(crossprod(gam_model_sim$family$data$R)))[1,2]
}

mean(unlist(corre_l))
quantile(unlist(corre_l), c(0.025, 0.975))

# Obtain confidence intervals for the SMOOTH model -------------------------------------------------
mod_smoo <- gam(list(log(CC_obs_loc1) ~ s(week_no, k = 50),
                     log(CC_obs_loc2) ~ s(week_no, k = 50)),
                family = mvn(d = 2), data = df_mod)

# Get Sigma
cov_mat_y = solve(crossprod(mod_smoo$family$data$R))
# Get the means 
mu_pred = predict(mod_smoo)
# Simulate 1
ysim = mvrnorm(n = 100, mu = mu_pred[1, ], Sigma = cov_mat_y)
# Simulate all 
ysim <- list()
for (i in 1:length(mu_pred[ ,1])) {
  ysim[[i]] = as.data.frame(mvrnorm(n = 100, mu = mu_pred[i, ], Sigma = cov_mat_y)) %>%
    rownames_to_column(var = "sim")
}
# Bind all simulations 
ysim_df <- bind_rows(ysim, .id = "week_no")

corre_l <- list()
for (sim_i in 1:length(unique(ysim_df$sim))) {
  # Add to df 
  df_mod_sim <- df_mod %>%
    mutate(loc1_sim = c(filter(ysim_df, sim == sim_i)$V1),
           loc2_sim = c(filter(ysim_df, sim == sim_i)$V2),
           loc1_sim = exp(loc1_sim),
           loc2_sim = exp(loc2_sim))
  
  # Run the model again 
  gam_model_sim <- gam(list(log(loc1_sim) ~ s(week_no, k = 50),
                            log(loc2_sim) ~ s(week_no, k = 50)),
                       family = mvn(d = 2), data = df_mod_sim)
  # Extract correlation 
  corre_l[[sim_i]] = cov2cor(solve(crossprod(gam_model_sim$family$data$R)))[1,2]
}

mean(unlist(corre_l))
quantile(unlist(corre_l), c(0.025, 0.975))

####################################################################################################
# END
####################################################################################################












