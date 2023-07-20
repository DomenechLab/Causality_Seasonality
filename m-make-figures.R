#######################################################################################################
# Make figures for all vignettes
#######################################################################################################

# Load packages -----------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))
save_plot <- F # Should all the plots be saved as a pdf? 

# Vignette on quasi-experiments -------------------------------------------

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
  theme(strip.background = element_blank()) + 
  labs(x = "Time (weeks)", y = "Total cases (per week per 100)")
print(pl_up)

# Lower panels: boxplots of climatic parameter estimates
tmp <- list(res_all[[1]]$pars_est, res_all[[2]]$pars_est)
names(tmp) <- loc_nm
tmp <- bind_rows(tmp, .id = "loc")

pl_low <- ggplot(data = tmp %>% filter(par %in% c("e_Te", "e_RH")), 
                 mapping = aes(x = par, y = mle)) + 
  geom_boxplot(outlier.color = "grey") + 
  geom_point(data = pars_true %>% filter(par %in% c("e_Te", "e_RH")), 
             mapping = aes(x = par, y = true), 
             color = "blue", 
             size = rel(3), shape = 17) +
  facet_wrap(~ loc, scales = "fixed", ncol = 2) + 
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "grey92"), 
        axis.line.y.left = element_line(color = "black"),
        strip.text = element_blank()) + 
  scale_x_discrete(labels = c("e_Te" = "Temperature", "e_RH" = "Relative humidity")) + 
  ylim(-0.5, 0.25) + 
  labs(x = "Climatic effect parameter", y = "Estimate")
print(pl_low)

# Assemble graph
pl_all <- pl_up / pl_low + 
  plot_annotation(tag_levels = "A")
print(pl_all)

ggsave(filename = "_figures/_main/vignette-quasi-experiments.pdf", 
       plot = pl_all, 
       width = 10, 
       height = 10)

# Make supp figures -------------------------------------------------------
tmp <- list(res_all[[1]]$pars_est, res_all[[2]]$pars_est)
names(tmp) <- loc_nm
tmp <- bind_rows(tmp, .id = "loc") %>% 
  filter(!(par %in% c("e_Te", "e_RH")))

# Transform parameters to natural scale
tmp$mle[tmp$par %in% c("R0", "alpha", "rho_k")] <- exp(tmp$mle[tmp$par %in% c("R0", "alpha", "rho_k")])
tmp$mle[tmp$par == "rho_mean"] <- plogis(tmp$mle[tmp$par == "rho_mean"]) 

pl_supp <- ggplot(data = tmp %>% filter(!(par %in% c("e_Te", "e_RH"))), 
                  mapping = aes(x = loc, y = mle)) + 
  geom_boxplot(outlier.color = "grey") + 
  geom_hline(data = pars_true %>% filter(!(par %in% c("e_Te", "e_RH"))),
             mapping = aes(yintercept = true),
             color = "blue", 
             linetype = "dashed") +
  facet_wrap(~ par, scales = "free", ncol = 2, 
             labeller = as_labeller(c("alpha" = "Waning rate (per week)", 
                                      "R0" = "Reproduction number", 
                                      "rho_k" = "Reporting overdispersion", 
                                      "rho_mean" = "Average reporting probability"))) + 
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "grey92"), 
        strip.background = element_blank()) + 
  scale_y_log10() + 
  labs(x = "Parameter", y = "Estimate")
print(pl_supp)

ggsave(filename = "_figures/_supp/vignette-quasi-experiments.pdf", 
       plot = pl_supp, 
       width = 10, 
       height = 8)

#######################################################################################################
# END
#######################################################################################################