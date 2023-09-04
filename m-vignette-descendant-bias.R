#######################################################################################################
# Run simulations for vignette on descendant bias
# Key point: many studies regress incidence on climatic variables, 
# But in fact climatic variables affect transmission, and the association between transmission and incidence can be complex
# Illustrate by running SIS and SIR models with climatic data in Madrid (weather station "LEVS") 
#######################################################################################################

# Load packages -----------------------------------------------------------
rm(list = ls())
source("f-Pred_RH.R")
Pred_RH <- Vectorize(FUN = Pred_RH)
source("f-CreateMod.R")
source("f-CreateClimData.R")
source("s-base_packages.R")
library("mgcv")
library("pomp")
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))
save_plot <- T # Should all the plots be saved as a pdf? 

# Set model parameters ----------------------------------------------------
parms <- c("mu" = 1 / 80 / 52, # Birth rate 
           "N" = 5e6, # Total population size
           "R0" = 1.25, # Basic reproduction no
           "sigma_beta" = 0, # SD of environmental noise (set to 0 for a deterministic process model)
           "e_Te" = -0.2, # Effect of Te on transmission 
           "e_RH" = -0.2, # Effect of RH on transmission
           "eps" = 1, # Fraction of infections conferring sterilizing immunity (1: SIR, 0: SIS)
           "alpha" = 0, # Rate of waning immunity
           "rho_mean" = 0.1, # Average reporting probability 
           "rho_k" = 0.04) # Reporting over-dispersion

# Load climatic data in a given location------------------------------------------------------
# Bogota (Colombia): weather station "SKBO", rho(Te, RH) = -0.08, CV(RH) = 0.07
# Barranquilla (Colombia): SKBQ, rho(Te, RH) = -0.09, CV(RH) = 0.06
# Gijon (Spain): LEAS, rho(Te, RH) = 0.29, CV(RH) = 0.08
# Rostock (Germany): Rostock, rho(Te, RH) = -0.69, CV(RH) = 0.085

e_Te <- parms["e_Te"]
e_RH <- parms["e_RH"]
loc_nm <- "Rostock"
if(save_plot) pdf(file = sprintf("_saved/_vignette_descendent_bias/vignette-descendant-bias-%s.pdf", loc_nm), 
                  width = 8, height = 8)

clim_dat <- CreateClimData(loc_nm = loc_nm, n_years = 10)

# Calculate seasonal term of transmission rate
clim_dat <- clim_dat %>% 
  mutate(beta_seas = exp(e_Te * (Te_norm - 1)) * exp(e_RH * (RH_pred_norm - 1))) %>% 
  select(loc, week_date, week_no, week, everything())

# Data in long format
clim_dat_long <- clim_dat %>% 
  pivot_longer(cols = -c(loc, week_date, week_no, week), names_to = "var", values_to = "value")

# Plot time series of climatic variables
pl <- ggplot(data = clim_dat_long %>% filter(var %in% c("Te", "Td", "RH", "RH_pred")), 
             mapping = aes(x = week_date, y = value)) + 
  geom_line() + 
  facet_wrap(~ var, scales = "free_y", ncol = 2) + 
  labs(x = "Time (weeks)", y = "Value", 
       title = sprintf("Climatic data, time plot (location: %s)", loc_nm))
print(pl)

# Plot measured and predicted RH
pl <- ggplot(data = clim_dat, mapping = aes(x = RH, y = RH_pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  labs(title = "Agreement between RH and RH_pred")
print(pl)

# Plot time series of measured and predicted RH
pl <- ggplot(data = clim_dat, mapping = aes(x = week_date, y = RH)) + 
  geom_line() + 
  geom_line(mapping = aes(x = week_date, y = RH_pred), color = "red") + 
  labs(title = "Agreement between RH and RH_pred (red curve)")
print(pl)

# Season plots of all climatic variables
pl <- ggplot(data = clim_dat_long %>% filter(var %in% c("Te_norm", "Td_norm", "RH_pred_norm")), 
             mapping = aes(x = isoweek(week_date), y = value, 
                           #color = factor(year(week_date)), 
                           group = factor(year(week_date))
             )) + 
  geom_line(color = "grey") + 
  geom_smooth(mapping = aes(x = isoweek(week_date), y = value, group = NULL, color = NULL)) + 
  facet_wrap(~ var, scales = "fixed", ncol = 2) + 
  labs(x = "Week no", y = "Rescaled value", title = "Climatic data, season plot")
print(pl)

# Plot seasonal forcing of transmission rate
pl <- ggplot(data = clim_dat, mapping = aes(x = week_date, y = beta_seas)) + 
  geom_line() + 
  labs(x = "Time (weeks)", title = "Seasonal forcing in transmission rate")
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
PompMod <- CreateMod(covars = covars, lin_bool_val = T)

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
sim <- simulate(object = PompMod, 
                nsim = 1, 
                params = p_mat, 
                format = "data.frame")

# sim <- trajectory(object = PompMod, 
#                   params = p_mat, 
#                   format = "data.frame")
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
    scale_x_continuous(breaks = 0:10)
  
  if(s == "S") {
    pl <- pl + labs(x = "Year", y = "Proportion", 
                    title = s, 
                    subtitle = sprintf("rho_mean=%.1f, rho_k=%.2f, N=%.1f M, eps=%.1f, e_Te=%.1f, e_RH=%.1f, s_beta=%.2f", 
                                       parms["rho_mean"], parms["rho_k"], parms["N"] / 1e6, parms["eps"], 
                                       parms["e_Te"], parms["e_RH"], parms["sigma_beta"]))
  } else {
    pl <- pl + labs(x = "Year", y = "Proportion", title = s)
  }
  
  print(pl)
}

# Run regressions ---------------------------------------------------------
n_reps <- 100L # No of stochastic replicates

sims_all <- simulate(object = PompMod, 
                     params = p_mat, 
                     nsim = n_reps, 
                     seed = 2186L, 
                     format = "data.frame")

sims_all <- sims_all %>% 
  separate(col = ".id", sep = "_", remove = T, into = c(".id", "rep")) %>% 
  rename("week_no" = "week") %>% 
  mutate(.id = as.integer(.id), 
         rep = as.integer(rep)) %>% 
  select(.id, rep, week_no, everything()) %>% 
  left_join(y = all_parms) %>% 
  left_join(y = clim_dat %>% select(week_no, matches("Te_norm|RH_pred_norm"), beta_seas)) %>% 
  arrange(.id, rep, week_no)

# Function to run GAM regression model
f_reg <- function(df) {
  
  # Run GAM model
  M_nb <- gam(formula = CC_obs ~ 1 + Te_norm_lag + RH_pred_norm_lag + s(week_no, k = 50), 
              family = nb(link = "log"), 
              data = arrange(df, week_no))
  
  # Extract regression coefficients
  M_est <- coef(M_nb)
  M_est_se <- summary(M_nb)$se
  R2 <- summary(M_nb)$r.sq
  
  # Correlation between smooth(time) and log(S*I)
  s_time <- predict(object = M_nb, type = "terms")
  s_time <- s_time[, "s(week_no)"]
  log_SI <- log(df$S * df$I)
  rho <- cor(x = s_time, y = log_SI, method = "pearson")
  
  # Return
  out <- data.frame(e_Te = M_est["Te_norm_lag"], 
                    e_Te_se = M_est_se["Te_norm_lag"], 
                    e_RH = M_est["RH_pred_norm_lag"], 
                    e_RH_se = M_est_se["RH_pred_norm_lag"],
                    rho_sTime_logSI = rho,
                    R2 = R2)
  
  return(out)
}

# Run regressions
sim_reg <- bake(file = sprintf("_saved/_vignette_descendent_bias/vignette-descendant-bias-regressions-%s.rds", loc_nm), 
                seed = 2186L, 
                expr = {
                  sims_all %>% 
                    group_nest(.id, rep, eps, R0, alpha) %>% 
                    mutate(reg = purrr::map(.x = data, .f = f_reg,  .progress = T)) %>% 
                    unnest(cols = "reg") %>% 
                    select(-data)
                })

sim_reg_long <- sim_reg %>% 
  pivot_longer(cols = e_Te:R2) %>% 
  mutate(R0_fac = R0 %>% factor() %>% fct_relabel(.fun = ~ paste0("R0 = ", .x)))

# Plots -------------------------------------------------------------------
# Plot point estimates; x-axis: waning rate, panels: R0
svars_plot <- c("e_Te", "e_RH", "R2", "rho_sTime_logSI")
for(s in svars_plot) {
  pl <- ggplot(data = sim_reg_long %>% filter(name == s), 
               mapping = aes(x = factor(1 / alpha / 52), y = value)) + 
    geom_boxplot() + 
    #facet_grid(R0 ~ factor(1 / alpha / 52), scales = "fixed") + 
    facet_wrap(~ R0_fac, scales = "fixed", ncol = 2) + 
    labs(x = "Average duration of immunity (years)", y = "Point estimate", title = s)
  
  if(s %in% c("e_Te", "e_RH")) {
    pl <- pl + 
      geom_hline(yintercept = coef(PompMod, s), color = "red", linetype = "dashed")
  } 
  
  print(pl)
}

# Plot point estimates; x-axis: R0, panels: waning rate
for(s in svars_plot) {
  pl <- ggplot(data = sim_reg_long %>% filter(name == s), 
               mapping = aes(x = factor(round(R0, 2)), y = value)) + 
    geom_boxplot() + 
    #facet_grid(R0 ~ factor(1 / alpha / 52), scales = "fixed") + 
    facet_wrap(~ factor(1 / alpha / 52), scales = "fixed", ncol = 2) + 
    labs(x = "R0", y = "Point estimate", title = s)
  
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
               mapping = aes(x = factor(1 / alpha / 52), y = value)) + 
    geom_boxplot() + 
    facet_wrap(~ R0_fac, scales = "fixed", ncol = 2) + 
    labs(x = "Average duration of immunity (years)", y = "Estimate SE", title = s)
  print(pl)
  
  pl <- ggplot(data = sim_reg_long %>% filter(name == s), 
               mapping = aes(x = factor(round(R0, 2)), y = value)) + 
    geom_boxplot() + 
    facet_wrap(~ factor(1 / alpha / 52), scales = "fixed", ncol = 2) + 
    labs(x = "R0", y = "Estimate SE", title = s)
  print(pl)
}

# Correlation between s(time) and R2
pl <- ggplot(data = sim_reg, 
             mapping = aes(x = R2, y = rho_sTime_logSI, 
                           color = abs(e_Te - parms["e_Te"]) + abs(e_RH - parms["e_RH"]))) + 
  geom_point(alpha = 0.5) + 
  scale_color_viridis(option = "rocket", direction = -1) + 
  theme(legend.position = "top") + 
  labs(x = "R2", y = "rho(sTime, log_SI)", color = "B(Te) + B(RH)", 
       title = "Association between smooth(time) and R2")
print(pl)

pl <- ggplot(data = sim_reg, 
             mapping = aes(x = R2, y = abs(e_Te - parms["e_Te"]) + abs(e_RH - parms["e_RH"]))) + 
  geom_point(color = "grey", alpha = 1) + 
  theme_classic() + 
  labs(x = "R2", y = "B(Te) + B(RH)", title = "Total absolute bias vs. R2")
print(pl)

pl <- ggplot(data = sim_reg, 
             mapping = aes(x = rho_sTime_logSI, y = abs(e_Te - parms["e_Te"]) + abs(e_RH - parms["e_RH"]))) + 
  geom_point(color = "grey", alpha = 1) + 
  theme_classic() + 
  labs(x = "rho(sTime, log_SI)", y = "B(Te) + B(RH)", title = "Total absolute bias vs. rho(sTime, log_SI)")
print(pl)

# Estimation performance --------------------------------------------------
est_perf <- sim_reg %>% 
  group_by(.id, eps, R0, alpha) %>% 
  summarise(e_Te_bias_mean = mean(abs(e_Te - parms["e_Te"])),
            e_Te_bias_sd = sd(abs(e_Te - parms["e_Te"])), 
            e_Te_se_mean = mean(e_Te_se), 
            e_Te_se_sd = sd(e_Te_se), 
            e_Te_pow = mean(e_Te + 1.96 * e_Te_se < 0), 
            e_RH_bias_mean = mean(abs(e_RH - parms["e_RH"])),
            e_RH_bias_sd = sd(abs(e_RH - parms["e_RH"])), 
            e_RH_se_mean = mean(e_RH_se), 
            e_RH_se_sd = sd(e_RH_se), 
            e_RH_pow = mean(e_RH + 1.96 * e_RH_se < 0), 
            R2_mean = mean(R2), 
            R2_sd = sd(R2)) %>% 
  ungroup() %>% 
  arrange(alpha, R0)

# Make main figures -------------------------------------------------------

# FIGURE 1: Time series of seasonal drivers and model simulations
# Plot of Te and RH
pl <- ggplot(data = clim_dat_long %>% filter(var %in% c("Te_norm", "RH_pred_norm")), 
             mapping = aes(x = isoweek(week_date), 
                           y = value, 
                           group = interaction(factor(year(week_date)), var),
                           color = var)) + 
  geom_line(alpha = 0.5) + 
  geom_smooth(mapping = aes(x = isoweek(week_date), y = value, group = var, color = var), se = F) + 
  scale_color_brewer(palette = "Set2", labels = c("RH", "Te")) + 
  theme_classic() + 
  theme(legend.position = c(0.1, 0.9)) + 
  labs(x = "Week no", y = "Renormalized value", color = "")
print(pl)

# Individual plots of Te, RH, and beta_seas
tmp <- clim_dat_long %>% 
  filter(var %in% c("Te_norm", "RH_pred_norm", "beta_seas"), week_no >= 1) %>% 
  mutate(var = factor(var, 
                      levels = c("Te_norm", "RH_pred_norm", "beta_seas"), 
                      labels = c("Temperature", "Relative humidity", "Transmission rate")))

pl_left <- ggplot(data = tmp, 
                  mapping = aes(x = week_no, 
                                y = value)) + 
  geom_line() + 
  geom_vline(xintercept = 52 * (0:10) + 1, color = "grey", alpha = 0.5) + 
  facet_wrap(~ var, ncol = 1, scales = "free_y") + 
  theme_classic() +
  theme(strip.background = element_blank(), 
        #strip.text = element_text(colour = "black", size = rel(1.2)), 
        panel.spacing = unit(1, "cm")) + 
  labs(x = "Time (weeks)", y = "Renormalized value")
print(pl_left)

# Plot of CC and S
tmp <- sim_long %>% 
  filter(1 / alpha == 52, state_var %in% c("S", "CC")) %>% 
  mutate(state_var = factor(state_var, levels = c("CC", "S"), 
                            labels = c("Incidence rate (% per week)", "Susceptible prevalence (%)")), 
         R0 = factor(R0))
levels(tmp$R0) <- paste0("R0 = ", levels(tmp$R0))


pl_right <- ggplot(data = tmp, 
                   mapping = aes(x = week_no, y = 1e2 * value / N, color = R0, linetype = state_var)) + 
  geom_line(linewidth = 0.5) + 
  facet_wrap(~factor(R0), scales = "free_y", ncol = 1) + 
  scale_y_sqrt() + 
  scale_color_viridis(discrete = T, option = "rocket", direction = -1, end = 0.75) +
  #scale_linetype_discrete(breaks = c("solid", "dotted")) + 
  #scale_y_log10() + 
  theme_classic() + 
  theme(strip.background = element_blank(), 
        legend.position = "top", 
        panel.spacing = unit(1, "cm")) + 
  guides(color = "none") + 
  labs(x = "Time (weeks)", y = "Value", color = expression(R[0]), linetype = "")
print(pl_right)

# Assemble plot

pl_all <- pl_left + pl_right + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 10))
print(pl_all)

ggsave(plot = pl_all, 
       filename = sprintf("_figures/_main/vignette-descendant-bias-main-%s.pdf", loc_nm), 
       width = 10, 
       height = 8)

# FIGURE 2: Parameter estimates from regression ---------------------------------------------------------------
tmp <- sim_reg %>% 
  filter(1 / alpha == 52) %>% 
  rename("Te" = "e_Te", "RH" = "e_RH") %>% 
  pivot_longer(cols = c("Te", "RH"), names_to = "var", values_to = "est")
  
pl <- ggplot(data = sim_reg %>% filter(1 / alpha == 52), 
             mapping = aes(x = e_Te, y = e_Te_se, color = R2, shape = factor(R0))) + 
  geom_vline(xintercept = e_Te, linetype = "dotted") + 
  geom_point() + 
  #geom_xsidedensity(mapping = aes(y = stat(density))) + 
  #facet_wrap(~ factor(R0)) + 
  scale_color_viridis(option = "turbo", direction = -1) +
  theme_classic() + 
  theme(legend.position = "top") + 
  labs(x = "Point estimate of temperature effect", y = "Standard error of estimate", 
       color = expression(R^2), shape = expression(R[0]))
print(pl)

pl2 <- ggplot(data = sim_reg %>% filter(1 / alpha == 52), 
             mapping = aes(x = e_RH, y = e_RH_se, color = R2, shape = factor(R0))) + 
  geom_vline(xintercept = e_RH, linetype = "dotted") + 
  geom_point() + 
  #facet_wrap(~ factor(R0)) + 
  scale_color_viridis(option = "turbo", direction = -1) +
  theme_classic() + 
  theme(legend.position = "top") + 
  labs(x = "Point estimate of relative humidity effect", y = "Standard error of estimate", 
       color = expression(R^2), shape = expression(R[0]))
print(pl2)

ggsave(plot = pl, 
       filename = sprintf("_figures/_main/vignette-descendant-bias-main-%s-Te.pdf", loc_nm), 
       width = 8, 
       height = 8)

ggsave(plot = pl2, 
       filename = sprintf("_figures/_main/vignette-descendant-bias-main-%s-RH.pdf", loc_nm), 
       width = 8, 
       height = 8)

# End statements ----------------------------------------------------------
if(save_plot) dev.off()

#######################################################################################################
# END
#######################################################################################################
