#######################################################################################################
# Run simulations for vignette on descendant bias
# Key point: many studies regress incidence on climatic variables, 
# But in fact climatic variables affect transmission, and the association between transmission and incidence can be complex
# Illustrate by running SIS and SIR models with climatic data in Madrid (weather station "LEVS") 
#######################################################################################################


# Load packages -----------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
source("f-Pred_RH.R")
Pred_RH <- Vectorize(FUN = Pred_RH)
source("f-CreateMod.R")
source("f-CreateClimData.R")
library(mgcv)
library(pomp)
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))

# Set model parameters ----------------------------------------------------
parms <- c("mu" = 1 / 80 / 52, # Birth rate 
           "N" = 5e6, # Total population size
           "R0" = 1.25, # Basic reproduction no
           "e_Te" = -0.2, # Effect of Te on transmission 
           "e_RH" = -0.2, # Effect of RH on transmission
           "eps" = 1, # Fraction of infections conferring sterilizing immunity (1: SIR, 0: SIS)
           "alpha" = 0, # Rate of waning immunity
           "rho_mean" = 0.1, # Average reporting probability 
           "rho_k" = 0.04) # Reporting over-dispersion

# Load climatic data in a given location------------------------------------------------------
# Bogota: weather station "SKBO"
# Madrid: weather station "LEVS"

e_Te <- parms["e_Te"]
e_RH <- parms["e_RH"]
loc_nm <- "SKBO"

clim_dat <- CreateClimData(loc_nm = loc_nm, n_years = 10)

# Calculate seasonal term of transmission rate
clim_dat <- clim_dat %>% 
  mutate(beta_seas = exp(e_Te * (Te_norm - 1)) * exp(e_RH * (RH_pred_norm - 1)))

# Data in long format
clim_dat_long <- clim_dat %>% 
  pivot_longer(cols = Te:RH_pred_norm, names_to = "var", values_to = "value")

# Plot time series of climatic variables
pl <- ggplot(data = clim_dat_long %>% filter(var %in% c("Te", "Td", "RH", "RH_pred")), 
             mapping = aes(x = week_date, y = value)) + 
  geom_line() + 
  facet_wrap(~ var, scales = "free_y", ncol = 2)
print(pl)

# Plot measured and predicted RH
pl <- ggplot(data = clim_dat, mapping = aes(x = RH, y = RH_pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")
print(pl)

# Plot time series of measured and predicted RH
pl <- ggplot(data = clim_dat, mapping = aes(x = week_date, y = RH)) + 
  geom_line() + 
  geom_line(mapping = aes(x = week_date, y = RH_pred), color = "red")
print(pl)

# Season plots of all climatic variables
pl <- ggplot(data = clim_dat_long %>% filter(var %in% c("Te_norm", "Td_norm", "RH_pred_norm")), 
             mapping = aes(x = isoweek(week_date), y = value, 
                           #color = factor(year(week_date)), 
                           group = factor(year(week_date))
             )) + 
  geom_line(color = "grey") + 
  geom_smooth(mapping = aes(x = isoweek(week_date), y = value, group = NULL, color = NULL)) + 
  facet_wrap(~ var, scales = "fixed", ncol = 2)
print(pl)

# Plot seasonal forcing of transmission rate
pl <- ggplot(data = clim_dat, mapping = aes(x = week_date, y = beta_seas)) + 
  geom_line()
print(pl)

# Add lagged variables for Te_norm and RH_norm
clim_dat <- clim_dat %>% 
  mutate(Te_norm_lag = lag(x = Te_norm, n = 1L, order_by = week_no), 
         RH_pred_norm_lag = lag(x = RH_pred_norm, n = 1L, order_by = week_no))

# Prepare covariate table
covars <- clim_dat %>% 
  select(week_date, week_no, Te, Td, RH_pred) %>% 
  arrange(week_date)

# Create pomp model -------------------------------------------------------
PompMod <- CreateMod(covars_df = covars, lin_bool_val = T)

# Set parameters
pomp::coef(PompMod, names(parms)) <- unname(parms)
base_pars <- pomp::coef(PompMod)

# Run simulations ---------------------------------------------------------
pomp::coef(PompMod, names(base_pars)) <- unname(base_pars)
rho_mean_val <- unname(coef(PompMod, "rho_mean"))
rho_k_val <- 0.04 # Reporting overdispersion
N_val <- unname(coef(PompMod, "N"))

R0_seq <- c(1.25, 2.5, 5)
alpha_seq <- c(0, 1 / 52, 1 / (2 * 52))
eps_seq <- 1

all_parms <- expand_grid(R0 = R0_seq, alpha = alpha_seq, eps = eps_seq) %>% 
  mutate(.id = seq_len(nrow(.)))

p_mat <- parmat(params = coef(PompMod), 
                nrep = nrow(all_parms))
p_mat["eps", ] <- all_parms$eps
p_mat["R0", ] <- all_parms$R0
p_mat["alpha", ] <- all_parms$alpha

# Simulations
sim <- trajectory(object = PompMod, 
                  params = p_mat, 
                  format = "data.frame")
sim <- sim %>% 
  mutate(.id = as.integer(.id)) %>% 
  left_join(y = all_parms) %>% 
  mutate(N_sim = S + I + R, 
         CC_obs = rnbinom(n = length(CC), mu = rho_mean_val * CC, size = 1 / rho_k_val)) %>% 
  rename(week_no = week) %>% 
  left_join(y = clim_dat %>% select(week_no, matches("Te_norm|RH_pred_norm"), beta_seas)) %>% 
  select(.id, eps, R0, alpha, week_no, S, I, R, N_sim, CC, CC_obs, everything())

sim_long <- sim %>% 
  pivot_longer(cols = S:CC_obs, names_to = "state_var", values_to = "value") %>% 
  mutate(state_var = factor(state_var, levels = c("S", "I", "R", "N_sim", "CC", "CC_obs")), 
         N = N_val)

# Plot simulations
svars_plot <- c("S", "N_sim", "CC", "CC_obs")

for(s in svars_plot) {
  pl <- ggplot(data = sim_long %>% filter(state_var == s), 
               mapping = aes(x = week_no / 52, y = value / N)) + 
    geom_line() + 
    #scale_y_sqrt() +
    facet_grid(R0 ~ factor(1 / alpha / 52), scales = "free_y") + 
    #facet_wrap(~ .id + state_var, scales = "free_y", ncol = nrow(all_parms), dir = "v") + 
    labs(x = "Year", y = "Proportion", title = s)
  print(pl)
}

# Run regressions ---------------------------------------------------------
n_reps <- 100L # No of stochastic replicates

# Function to run regressions
f_reg <- function(df, n = 5L) {
  
  df <- df %>% 
    arrange(week_no) 
  
  M_nb <- vector(mode = "list", length = n)
  
  # Create observations and fit models
  for(i in 1:n) {
    df[[paste0("CC_obs_", i)]] <- rnbinom(n = length(df$CC), mu = rho_mean_val * df$CC, size = 1 / rho_k_val)
    df[[paste0("CC_obs_lag_", i)]] <-  lag(x = df[[paste0("CC_obs_", i)]], n = 1L, order_by = df$week_no)
    
    df_cur <- df %>% select("week_no", Te_norm_lag, RH_pred_norm_lag, paste0(c("CC_obs_", "CC_obs_lag_"), i))
    colnames(df_cur) <- c("week_no", "Te", "RH", "CC_obs", "CC_obs_lag")
    M_nb[[i]] <- gam(formula = CC_obs ~ 1 + log(CC_obs_lag + 1) + Te + RH + s(week_no, k = 10), 
                     family = nb(link = "log"), 
                     data = df_cur)
  }
  
  # Extract regression coefficients
  M_est <- sapply(M_nb, coef)
  M_est_se <- sapply(M_nb, function(y) summary(y)$se)
  R2 <- sapply(M_nb, function(y) summary(y)$r.sq)
  out <- data.frame(e_Te = M_est["Te", ], 
                    e_Te_se = M_est_se["Te", ], 
                    e_RH = M_est["RH", ], 
                    e_RH_se = M_est_se["RH", ], 
                    R2 = R2)
  
  return(out)
}

# Run regressions
sim_reg <- bake(file = sprintf("_saved/vignette-descendant-bias-regressions-%s.rds", loc_nm), 
                seed = 2186L, 
                expr = {
                  sim %>% 
                    group_nest(.id, eps, R0, alpha) %>% 
                    mutate(reg = purrr::map(.x = data, .f = f_reg, n = n_reps,  .progress = T)) %>% 
                    unnest(cols = "reg") %>% 
                    select(-data)
                })
 
sim_reg_long <- sim_reg %>% 
  pivot_longer(cols = e_Te:R2) %>% 
  mutate(R0_fac = R0 %>% factor() %>% fct_relabel(.fun = ~ paste0("R0 = ", .x)))



# Plot point estimates
svars_plot <- c("e_Te", "e_RH", "R2")
for(s in svars_plot) {
  pl <- ggplot(data = sim_reg_long %>% filter(name == s), 
               mapping = aes(x = factor(round(alpha, 2)), y = value)) + 
    geom_boxplot() + 
    #facet_grid(R0 ~ factor(1 / alpha / 52), scales = "fixed") + 
    facet_wrap(~ R0_fac, scales = "fixed", ncol = 2) + 
    labs(x = "Waning rate (per year)", y = "Point estimate", title = s)
  
  if(s %in% c("e_Te", "e_RH")) {
    pl <- pl + 
      geom_hline(yintercept = coef(PompMod, s), color = "red", linetype = "dashed")
  } 
  
  print(pl)
}

# Plot precision of estimates
svars_plot <- c("e_Te_se", "e_RH_se")
for(s in svars_plot) {
  pl <- ggplot(data = sim_reg_long %>% filter(name == s), 
               mapping = aes(x = factor(round(alpha, 2)), y = value)) + 
    geom_boxplot() + 
    facet_wrap(~ R0_fac, scales = "fixed", ncol = 2) + 
    labs(x = "Waning rate (per year)", y = "Estimate SE", title = s)
  
  print(pl)
}

#######################################################################################################
# END
#######################################################################################################
