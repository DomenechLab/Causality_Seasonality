#######################################################################################################
# Make figures for all vignettes
#######################################################################################################

# Load packages -----------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
source("f-Pred_RH.R")
Pred_RH <- Vectorize(FUN = Pred_RH)
source("f-CreateClimData.R")
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))
save_plot <- F # Should all the plots be saved as a pdf? 

#######################################################################################################
# VIGNETTE ON QUASI-EXPERIMENTS
#######################################################################################################

# Load climatic data ------------------------------------------------------
loc_nm <- c("Bogota", "Rostock")
clim_dat <- vector(mode = "list", length = 2)
names(clim_dat) <- loc_nm

clim_dat[["Bogota"]] <- CreateClimData(loc_nm = "SKBO", n_years = 10)
clim_dat[["Rostock"]] <- CreateClimData(loc_nm = "Rostock", n_years = 10)

clim_dat <- clim_dat %>% 
  bind_rows(.id = "loc") %>% 
  pivot_longer(cols = -c("loc", "week_date", "week_no", "week"), 
               names_to = "clim_var", 
               values_to = "value")

pl <- ggplot(data = clim_dat %>% filter(clim_var %in% c("Te_norm", "RH_pred_norm")), 
             mapping = aes(x = week, 
                           y = value, 
                           group = interaction(factor(year(week_date)), clim_var),
                           color = clim_var)) + 
  geom_line(alpha = 0.5) + 
  geom_smooth(mapping = aes(x = week, y = value, group = clim_var, color = clim_var), se = F) + 
  facet_wrap(~ loc, scales = "free_y", ncol = 2) + 
  scale_color_brewer(palette = "Set2", labels = c("RH", "Te")) + 
  theme_classic() + 
  theme(legend.position = "top", 
        strip.background = element_blank(), 
        strip.text = element_text(colour = "black", size = rel(1))) + 
  labs(x = "Week no", y = "Renormalized value", color = "")
print(pl)

ggsave(plot = pl, 
       filename = "_figures/_supp/vignette-quasi-experiments-s1.pdf", 
       width = 10, 
       height = 8)

# Load simulations and estimations  ---------------------------------------
loc_nm <- c("Rostock", "Bogota") # Names of locations
res_all <- vector(mode = "list", length = length(loc_nm)) # List to contain all the results
names(res_all) <- loc_nm
pl1 <- pl2 <- vector(mode = 'list', length = length(loc_nm))
names(pl1) <- names(pl2) <- loc_nm

# Load results
res_all[["Rostock"]] <- readRDS(file = "_saved/vignette-quasi-experiments-Rostock-all.rds")
res_all[["Bogota"]] <- readRDS(file = "_saved/vignette-quasi-experiments-SKBO-all.rds")
N_val <- res_all[[1]]$sim_det$N_sim[1] # Population size

# True parameter values
pars_true <- res_all[[1]]$pars_est %>% 
  filter(sim_no == 1) %>% 
  select(par, true)

# Convert to natural scale
pars_true$true[pars_true$par %in% c("R0", "alpha", "rho_k")] <- exp(pars_true$true[pars_true$par %in% c("R0", "alpha", "rho_k")])
pars_true$true[pars_true$par == "rho_mean"] <- plogis(pars_true$true[pars_true$par == "rho_mean"])

# Convert stochastic simulations to long format
for(i in seq_along(loc_nm)) {
  res_all[[i]]$sim_stoch <- res_all[[i]]$sim_stoch %>% 
    pivot_longer(-c(".id", "week"), names_to = "var", values_to = "value")
}


# Make main figure --------------------------------------------------------
# Upper panels: time series of CC_t, for deterministic and 10 random stochastic simulations
id_sims <- pomp::freeze(expr = sample(x = 1:max(res_all[[1]]$sim_stoch$.id), size = 10, replace = F), 
                        seed = 2186L) # random indices for 10 stochastic simulations

tmp <- list(res_all[[1]]$sim_stoch, res_all[[2]]$sim_stoch)
tmp2 <- list(res_all[[1]]$sim_det, res_all[[2]]$sim_det)
names(tmp) <- names(tmp2) <- loc_nm
tmp <- bind_rows(tmp, .id = "loc")
tmp2 <- bind_rows(tmp2, .id = "loc")

pl_up <- ggplot(data = tmp %>% filter(.id %in% id_sims, var == "CC"), 
                mapping = aes(x = week, y = 1e2 * value / N_val, group = .id)) + 
  geom_line(color = "grey", alpha = 0.5) + 
  geom_line(data = tmp2, mapping = aes(x = week, y = 1e2 * CC / N_val), color = "black") +
  facet_wrap(~ loc, scales = "free_y", ncol = 2) + 
  theme_classic() + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = rel(1), colour = "black")) + 
  labs(x = "Time (weeks)", y = "Total cases (per week per 100)")
print(pl_up)

# Lower panels: boxplots of climatic parameter estimates
tmp <- list(res_all[[1]]$pars_est, res_all[[2]]$pars_est)
names(tmp) <- loc_nm
tmp <- bind_rows(tmp, .id = "loc")
tmp <- tmp %>% 
  filter(par %in% c("e_Te", "e_RH")) %>% 
  mutate(par = factor(par, levels = c("e_Te", "e_RH"), labels = c("Temperature", "Relative humidity")))

pl_low <- ggplot(data = tmp, 
                 mapping = aes(x = mle, color = loc, fill = loc)) + 
  geom_vline(xintercept = pars_true$true[pars_true$par == "e_Te"], linetype = "dotted") + 
  geom_density(alpha = 0.1, adjust = 2) + 
  geom_rug(alpha = 0.5) + 
  facet_wrap(~ par, scales = "fixed", ncol = 2) + 
  scale_color_brewer(palette = "Set2") + 
  scale_fill_brewer(palette = "Set2") + 
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = rel(1), colour = "black"),
        legend.position = c(0.5, 0.8)) + 
  labs(x = "Parameter estimate", y = "Density", color = "", fill = "")
print(pl_low)

# Assemble graph
pl_all <- pl_up / pl_low + 
  plot_annotation(tag_levels = "A")
print(pl_all)

ggsave(filename = "_figures/_main/vignette-quasi-experiments-m1.pdf", 
       plot = pl_all, 
       width = 10, 
       height = 8)

# Make supp figures -------------------------------------------------------
tmp <- list(res_all[[1]]$pars_est, res_all[[2]]$pars_est)
names(tmp) <- loc_nm
tmp <- bind_rows(tmp, .id = "loc") %>% 
  filter(!(par %in% c("e_Te", "e_RH")))

# Transform parameters to natural scale
tmp$mle[tmp$par %in% c("R0", "alpha", "rho_k")] <- exp(tmp$mle[tmp$par %in% c("R0", "alpha", "rho_k")])
tmp$mle[tmp$par == "rho_mean"] <- plogis(tmp$mle[tmp$par == "rho_mean"]) 

pl_supp <- ggplot(data = tmp, 
                  mapping = aes(x = mle, color = loc, fill = loc)) + 
  geom_vline(data = pars_true %>% filter(!(par %in% c("e_Te", "e_RH"))),
             mapping = aes(xintercept = true),
             color = "black", 
             linetype = "dotted") +
  geom_density(alpha = 0.1, adjust = 2) + 
  geom_rug(alpha = 0.5) + 
  facet_wrap(~ par, 
             scales = "free", 
             ncol = 2, 
             labeller = as_labeller(c("alpha" = "Waning rate (per week)", 
                                      "R0" = "Reproduction number", 
                                      "rho_k" = "Reporting overdispersion", 
                                      "rho_mean" = "Average reporting probability"))) + 
  scale_color_brewer(palette = "Set2") + 
  scale_fill_brewer(palette = "Set2") + 
  theme_classic() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = rel(1), colour = "black"), 
        legend.position = "top") + 
  scale_x_log10() + 
  labs(x = "Parameter estimate", y = "Density", color = "", fill = "")
print(pl_supp)

ggsave(filename = "_figures/_supp/vignette-quasi-experiments-s2.pdf", 
       plot = pl_supp, 
       width = 10, 
       height = 8)

#######################################################################################################
# END
#######################################################################################################