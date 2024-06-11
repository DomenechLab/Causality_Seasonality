####################################################################################################
# Run simulations to check the model with coupling
####################################################################################################

# Load packages ------------------------------------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
source("f-CreateMod2coup.R")
library(pomp)
library(glue)
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))

# Load climatic data in a given location------------------------------------------------------------
# Choose between Spain and Colombia
coun_name <- "Colombia"
# Choose 2 locations
if(coun_name == "Spain") covarsnam_2loc <- c("LEVS", "LEAS") # Select only Madrid and Gijon 
if(coun_name == "Colombia") covarsnam_2loc <- c("SKBO", "SKPS") # Select only Bogota and Pasto 

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

# Create pomp model --------------------------------------------------------------------------------
PompMod <- CreateMod(covars = covars, lin_bool_val = T)

# Consider different R0 and p (coupling) values with and without climate forcing
all_parms <- expand_grid(R0 = c(2.5), 
                         p = c(0, 0.001, 0.01, 0.1, 0.5, 0.75, 1), 
                         e_Te = c(0, -0.2),
                         sigma_beta = c(0, 0.02, 0.2, 0.3, 0.4, 0.5, 1)) %>% 
  mutate(e_RH = e_Te) %>%
  mutate(.id = seq_len(nrow(.)))

# Create list of all parameters
p_mat <- parmat(params = coef(PompMod), nrep = nrow(all_parms))
p_mat["R0", ] <- all_parms$R0
p_mat["p", ] <- all_parms$p
p_mat["e_Te", ] <- all_parms$e_Te
p_mat["e_RH", ] <- all_parms$e_RH
p_mat["sigma_beta", ] <- all_parms$sigma_beta

# Simulate 
sim <- simulate(object = PompMod, params = p_mat, format = "data.frame") %>%
  mutate(.id = as.integer(.id)) %>% 
  left_join(y = all_parms) %>%
  mutate(N1 = S1 + I1 + R1,
         N2 = S2 + I2 + R2) %>%
  relocate(N1, N2) %>%
  pivot_longer(cols = c(S1, I1, R1, N1, CC1, CC_obs1, S2, I2, R2, CC2, N2, CC_obs2))

# Plot trajectories (no extra-demographic stochasticity)
sim %>%  
  filter(sigma_beta == 0, e_Te == 0) %>%
  separate(name, into = c("name", "Patch"), sep = -1) %>%
  filter(name != "CC") %>%
  mutate(name = factor(name, levels = c("S", "I", "R", "N", "CC_obs"))) %>%
  ggplot(aes(x = week, y = value, color = Patch, linetype = Patch)) + 
    geom_line() + 
    scale_color_viridis_d(end = 0.75) +
    scale_y_continuous("Number") +
    scale_x_continuous("Week") +
    facet_grid(name~p, scales = "free") + 
    ggtitle("Without climate forcing and no extra-demographic stochasticity") + 
sim %>%  
  filter(sigma_beta == 0, e_Te == -0.2) %>%
  separate(name, into = c("name", "Patch"), sep = -1) %>%
  filter(name != "CC") %>%
  mutate(name = factor(name, levels = c("S", "I", "R", "N", "CC_obs"))) %>%
  ggplot(aes(x = week, y = value, color = Patch, linetype = Patch)) + 
  geom_line() + 
  scale_color_viridis_d(end = 0.75) +
  scale_y_continuous("Number") +
  scale_x_continuous("Week") +
  facet_grid(name~p, scales = "free") + 
  ggtitle("With climate forcing and no extra-demographic stochasticity") + 
plot_layout(ncol = 1)

# Plot simulations (with extra-demographic stochasticity)
sim %>%  
  filter(sigma_beta == 0.02, e_Te == 0) %>%
  separate(name, into = c("name", "Patch"), sep = -1) %>%
  filter(name != "CC") %>%
  mutate(name = factor(name, levels = c("S", "I", "R", "N", "CC_obs"))) %>%
  ggplot(aes(x = week, y = value, color = Patch, linetype = Patch)) + 
  geom_line() + 
  scale_color_viridis_d(end = 0.75) +
  scale_y_continuous("Number") +
  scale_x_continuous("Week") +
  facet_grid(name~p, scales = "free") + 
  ggtitle("Without climate forcing and no extra-demographic stochasticity") + 
sim %>%  
  filter(sigma_beta == 0.02, e_Te == -0.2) %>%
  separate(name, into = c("name", "Patch"), sep = -1) %>%
  filter(name != "CC") %>%
  mutate(name = factor(name, levels = c("S", "I", "R", "N", "CC_obs"))) %>%
  ggplot(aes(x = week, y = value, color = Patch, linetype = Patch)) + 
  geom_line() + 
  scale_color_viridis_d(end = 0.75) +
  scale_y_continuous("Number") +
  scale_x_continuous("Week") +
  facet_grid(name~p, scales = "free") + 
  ggtitle("With climate forcing and no extra-demographic stochasticity") + 
plot_layout(ncol = 1)

# Check for sigma beta 
sim %>%  
  filter(p == 0) %>%
  separate(name, into = c("name", "Patch"), sep = -1) %>%
  filter(name == "CC_obs") %>%
  mutate(name = factor(name, levels = c("S", "I", "R", "N", "CC_obs"))) %>%
  ggplot(aes(x = week, y = value, color = Patch, linetype = Patch)) + 
  geom_line() + 
  scale_color_viridis_d(end = 0.75) +
  scale_y_continuous("Number") +
  scale_x_continuous("Week") +
  facet_wrap(sigma_beta~e_Te, scales = "free", ncol = 2) + 
  ggtitle("With and without climate forcing, no coupling, different values for sigma_beta") + 
  plot_layout(ncol = 1)

# Check 
sim %>%  
  filter(e_Te == 0) %>%
  separate(name, into = c("name", "Patch"), sep = -1) %>%
  filter(name == "CC_obs") %>%
  mutate(name = factor(name, levels = c("S", "I", "R", "N", "CC_obs"))) %>%
  ggplot(aes(x = week, y = value, color = Patch, linetype = Patch)) + 
  geom_line() + 
  scale_color_viridis_d(end = 0.75) +
  scale_y_continuous("Number") +
  scale_x_continuous("Week") +
  facet_wrap(p~sigma_beta, scales = "free", ncol = 7) + 
  ggtitle("Without climate forcing, different values for p and sigma_beta") + 
  plot_layout(ncol = 1)

#######################################################################################################
# END
#######################################################################################################
