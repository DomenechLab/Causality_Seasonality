
# Load packages -----------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
theme_set(theme_classic())

# MAPS =============================================================================================
# Map of globe

world <- map_data("world")
colombia <- map_data("world") %>%
  filter(region == "Colombia")

globe_colombia <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill = "grey40") +
  geom_polygon(data = colombia, aes(x=long, y = lat, group = group), 
               fill = "firebrick2") + 
  geom_hline(aes(yintercept = 0), size = 0.5) + 
  coord_map("ortho", orientation = c(0, -120, 0)) + 
  theme_minimal() + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())

germany <- map_data("world") %>%
  filter(region == "Germany")

globe_germany <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill = "grey40") +
  geom_polygon(data = germany, aes(x=long, y = lat, group = group), 
               fill = "firebrick2") + 
  geom_hline(aes(yintercept = 0), size = 0.5) + 
  coord_map("ortho", orientation = c(0, -40, 0)) + 
  theme_minimal() + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())

globe_germany + globe_colombia

# Figure 4 =========================================================================================

sim_long <- readRDS(glue("_saved/_vignette_total_effect/sim_te_effect_Pasto.rds")) %>%
  bind_rows(readRDS(glue("_saved/_vignette_total_effect/sim_te_effect_Lübeck.rds")), 
            .id = "loc_tidy_nm") %>%
  mutate(loc_tidy_nm = case_when(loc_tidy_nm == 1 ~ "Pasto",
                                 loc_tidy_nm == 2 ~ "Lübeck"))

col <- c("#A7B9CA","#134471")

fun_main_pl <- function(loc = "Pasto", legend = T) {
  pl1a <- sim_long %>%
    filter(state_var == "CC") %>%
    filter(week_no <= 312) %>%
    filter(R0 == 1.25) %>% 
    filter(loc_tidy_nm == loc) %>%
    ggplot(aes(x = week_no, y = beta_seas, color = Te_effect, group = Te_effect)) + 
    geom_line(size = 0.8) + 
    scale_x_continuous("Time (weeks)") + 
    scale_y_continuous("Renormalized \n transmission rate") + 
    scale_color_manual("Effect of temperature", values = col,
                       labels = c(expression(paste("Indirect effect (", delta[Te], " = 0, ", delta[RH], " = -0.2)")),
                                  expression(paste("Total effect (", delta[Te], " = -0.2, ", delta[RH], " = -0.2)")))) +
    facet_wrap(.~loc_tidy_nm, scales = "free", ncol = 2) + 
    theme(strip.text = element_text(size = rel(1), colour = "black"),
          legend.position = "top",
          legend.direction = "vertical",
          legend.justification = "right",
          legend.title.align = 1,
          legend.text.align = 1,
          strip.background = element_blank(),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          legend.background = element_rect(fill='transparent')
          )
  
  if(legend == F) {pl1a <- pl1a + theme(legend.position = "none")}
  
  pl1b <- sim_long %>%
    filter(state_var == "CC") %>%
    filter(R0 == 1.25) %>% 
    mutate(Te_effect = case_when(str_detect(Te_effect, "Indirect") ~ "Indirect effect",
                                 str_detect(Te_effect, "Total") ~ "Total effect")) %>%
    filter(loc_tidy_nm == loc) %>%
    ggplot(aes(x = Te, y = beta_seas, color = Te_effect, group = Te_effect)) + 
    geom_point() + 
    scale_x_continuous("Temperature (°C)") + 
    scale_y_continuous("Renormalized \n transmission rate") + 
    facet_wrap(Te_effect~., scales = "free", ncol = 4) + 
    scale_color_manual("Effect of temperature", values = col) + 
    theme(strip.text = element_text(size = rel(1), colour = "black"),
          legend.position = "none",
          strip.background = element_blank(),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          legend.background = element_rect(fill='transparent')
          )
  
  pl1 <- pl1a / pl1b
  return(pl1) 
}

pl1 <- fun_main_pl(loc = "Lübeck", legend = F) | fun_main_pl(loc = "Pasto") +   
  plot_annotation(theme = theme(plot.background = element_rect(color  = 'transparent', fill ="transparent")))

pl1

# ggsave(pl1, file = glue("_figures/00_Fig_total-effect.pdf"),
#        height = 7, width = 9)


# Figure 3 =========================================================================================

sims_stoch <- readRDS(glue("_presentationfigures/vignette-quasi-experiments-sims_stoch.rds"))
sims_det <- readRDS(glue("_presentationfigures/vignette-quasi-experiments-sims_det.rds"))
pars_main <- readRDS(glue("_presentationfigures/vignette-quasi-experiments-pairs_main.rds"))
id_sims <- readRDS(glue("_presentationfigures/vignette-quasi-experiments-id_sims.rds"))

pl_up_bogota <- ggplot(data = sims_stoch %>% filter(.id %in% id_sims, var == "CC", loc == "Bogota"), 
                mapping = aes(x = week, y = 1e2 * value / 5e+06, group = .id)) + 
  geom_line(color = "grey", alpha = 0.5) + 
  geom_line(data = filter(sims_det, loc == "Bogota"), mapping = aes(x = week, y = 1e2 * CC / 5e+06), color = "black",
            size = 0.8) +
  facet_wrap(~ loc, scales = "free_y", nrow = 1) + 
  theme_classic() + 
  theme(strip.text = element_text(size = rel(1), colour = "black"),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent')) + 
  labs(x = "Time (weeks)", y = "Total cases \n (per week per 100)")
print(pl_up_bogota)

pl_up_lubeck <- ggplot(data = sims_stoch %>% filter(.id %in% id_sims, var == "CC", loc == "Lübeck"), 
                       mapping = aes(x = week, y = 1e2 * value / 5e+06, group = .id)) + 
  geom_line(color = "grey", alpha = 0.5) + 
  geom_line(data = filter(sims_det, loc == "Lübeck"), mapping = aes(x = week, y = 1e2 * CC / 5e+06), color = "black",
            size = 0.8) +
  facet_wrap(~ loc, scales = "free_y", nrow = 1) + 
  theme_classic() + 
  theme(strip.text = element_text(size = rel(1), colour = "black"),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent')) + 
  labs(x = "Time (weeks)", y = "Total cases \n (per week per 100)")
print(pl_up_lubeck)

pl_low_bogota <- ggplot(data = filter(pars_main, loc == "Bogota"), 
                 mapping = aes(x = mle, color = par, fill = par, alpha = par)) + 
  geom_density(adjust = 2, color = "white") + 
  geom_rug(alpha = 0.7) + 
  geom_vline(xintercept = -0.2, linetype = "dotted") + 
  facet_wrap(. ~ par, scales = "free_y", ncol = 2) + 
  scale_color_manual(values = c("#F2B336", "grey20")) + 
  scale_fill_manual(values = c("#F2B336", "grey20")) + 
  scale_alpha_manual(values = c(0.9, 0.4)) + 
  scale_x_continuous("Parameter estimate", limits = c(-0.5, 0.75)) + 
  scale_y_continuous("Density", limits = c(0, 10)) + 
  guides(alpha = "none") + 
  theme_classic() +
  theme(strip.text = element_text(size = rel(1), colour = "black"),
        legend.position = "none",
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'))
print(pl_low_bogota)

pl_low_lubeck <- ggplot(data = filter(pars_main, loc == "Lübeck"), 
                        mapping = aes(x = mle, color = par, fill = par, alpha = par)) + 
  geom_density(adjust = 2, color = "white") + 
  geom_rug(alpha = 0.7) + 
  geom_vline(xintercept = -0.2, linetype = "dotted") + 
  facet_wrap(. ~ par, scales = "free_y", ncol = 2) + 
  scale_color_manual(values = c("#F2B336", "grey20")) + 
  scale_fill_manual(values = c("#F2B336", "grey20")) + 
  scale_alpha_manual(values = c(0.9, 0.4)) + 
  scale_x_continuous("Parameter estimate", limits = c(-0.5, 0.75)) + 
  scale_y_continuous("Density", limits = c(0, 28)) + 
  guides(alpha = "none") + 
  theme_classic() +
  theme(strip.text = element_text(size = rel(1), colour = "black"),
        legend.position = "none",
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'))
print(pl_low_lubeck)

# Assemble graph
(pl_up_bogota / pl_low_bogota) | fun_main_pl(loc = "Pasto", legend = F) 

(pl_up_lubeck / pl_low_lubeck) | fun_main_pl(loc = "Lübeck", legend = F) 

#10*5

pl_all <- pl_up / pl_low
print(pl_all)


















# Figure 1 =========================================================================================

# Panel A: Individual plots of Te, RH, and beta_seas
tmpa <- clim_dat_long %>% 
  filter(var %in% c("Te_norm", "RH_pred_norm", "beta_seas"), week_no >= 1) %>% 
  mutate(var = factor(var, 
                      levels = c("Te_norm", "RH_pred_norm", "beta_seas"), 
                      labels = c("Temperature", "Relative humidity", "Transmission rate"))) %>%
  filter(week_no <= 312)

pl_A <- ggplot(data = tmpa, 
               mapping = aes(x = week_no, 
                             y = value)) + 
  geom_line(size = 0.8) + 
  geom_vline(xintercept = 52 * (0:6) + 1, color = "grey", alpha = 0.5) + 
  facet_wrap(~ var, ncol = 1, scales = "free_y") + 
  theme_classic() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent')) + 
  labs(x = "Time (weeks)", y = "Renormalized value")
print(pl_A)

# PANEL B: Plot of CC and S
tmpb <- sim_long %>% 
  filter(1 / alpha == 52, state_var %in% c("S", "CC")) %>% 
  mutate(state_var = factor(state_var, levels = c("S", "CC"), 
                            labels = c("Susceptible prevalence", "Incidence rate")), 
         R0 = factor(R0)) %>%
  filter(week_no <= 312)
levels(tmpb$R0) <- paste0("R0 = ", levels(tmpb$R0))


pl_B <- ggplot(data = tmpb, 
               mapping = aes(x = week_no, y = 1e2 * value / N, linetype = state_var)) + 
  geom_line(size = 0.6) + 
  facet_wrap(~factor(R0), scales = "free_y", ncol = 1) + 
  scale_y_sqrt() + 
  scale_linetype_manual(values = c(11, "solid")) + 
  theme_classic() + 
  theme(legend.position = "bottom", 
        #legend.direction = "vertical",
        legend.text = element_text(size = 10),
        panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent')) + 
  guides(color = "none") + 
  labs(x = "Time (weeks)", y = "Value", color = expression(R[0]), linetype = "")
print(pl_B)

tmpc <- sim_reg_all %>% 
  filter(1 / alpha == 52, model == "smooth") %>% 
  mutate(R0 = factor(R0))
levels(tmpc$R0) <- paste0("R0 = ", levels(tmpc$R0))

# PANEL C: distribution of estimates
pl_C <- ggplot(data = tmpc, 
               mapping = aes(x = e_Te, color = e_Te_se)) + 
  geom_vline(xintercept = -0.2, linetype = "dotted", color = "black") + 
  geom_density(linewidth = 0.8, color = "black") + 
  geom_rug() + 
  facet_wrap(~ R0, ncol = 1) + 
  scale_colour_gradient2(low = "#FDF5E6", mid = "#F2B336", high = "#1E1607",
                         midpoint = 0.08) + 
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0)) + 
  theme_classic() + 
  theme(legend.position = "right", 
        panel.spacing = unit(1, "lines"), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent')) + 
  labs(x = "Point estimate of \n temperature effect",
       y = "Density",
       color = "Standard \nerror")
print(pl_C)

# Assemble plots ----------------------------------------------------------
pl_all <- pl_A + pl_B + pl_C +  
  plot_layout(widths = c(1,1,0.5))
print(pl_all)


saveRDS(pl_all, "pl1.rds")

ggsave(plot = pl_all, 
         filename = "presentation_pl1.svg", 
         width = 10, 
         height = 6)

# Figure 2 =========================================================================================

sncf_df <- bind_rows(lapply(simulations_l, filter, sim_id == 1)) %>%
  filter(state_var == "CC_obs") %>%
  filter(.id == 1) %>%
  mutate(inc = value/N) %>% 
  left_join(spat_dat[c("loc","lon","lat","loc_city_name")]) %>%
  select(loc, lon, lat, week_no, inc) %>%
  arrange(week_no) %>%
  pivot_wider(names_from = week_no,
              values_from = inc)

# Estimate the nonparametric (cross-)correlation function
sncf_result0 <- Sncf(x = sncf_df$lon, y = sncf_df$lat, z = sncf_df[,4:523], resamp=1000, latlon = TRUE)

# Get labels for the plots
all_params_labs <- all_parms %>%
  mutate(mod = case_when(eps == 0 ~ "SIS",
                         eps == 1 ~ "SIR"),
         alpha = 1/(52*alpha)) %>%
  filter(.id == 1)

fun_sncf_pl <- function(sncf_result = sncf_result0) {data.frame(x = c(sncf_result$real$predicted$x),
                                                                y = c(sncf_result$real$predicted$y),
                                                                ylow = c(sncf_result$boot$boot.summary$predicted$y["0.025", ]),
                                                                yhig = c(sncf_result$boot$boot.summary$predicted$y["0.975", ])) %>%
    mutate(ymean = mean(y)) %>%
    mutate(yhig = case_when(yhig > 1 ~ 1,
                            .default = yhig),
           ylow = case_when(ylow < -0.1 ~ -0.1,
                            .default = ylow))}

sncf_pl0 <- fun_sncf_pl(sncf_result = sncf_result0)

# Plot simulations
sim_pl <- bind_rows(lapply(simulations_l, filter, sim_id == 1)) %>%
  filter(week_no > 200 & week_no <= 500) %>%
  mutate(week_no = week_no - 200) %>%
  filter(state_var == "CC_obs") %>%
  filter(.id == 1) %>%
  left_join(spat_dat[c("loc","lon","lat","loc_city_name")]) %>%
  mutate(loc_lab = glue("{loc_city_name} ({loc})")) %>%
  group_by(.id, loc) %>%
  mutate(inc = value/N,
         inc_minmax = (inc - min(inc))/ (max(inc) - min(inc))) %>%
  ggplot(aes(fill = inc_minmax, x = week_no, y = reorder(loc_lab, lat))) + 
  geom_tile(linewidth = 0.0001) +
  # scale_fill_viridis_c("Scaled incidence",
  #                      labels = function(x) format(round(x, 1), nsmall = 1)) + 
  #scale_fill_gradient2(low = "firebrick2", mid = "#F9F5D2", high = "#134471", midpoint = 0.50) +
   scale_fill_gradient("Scaled incidence", low = "white", high = "#134471",
                       labels = function(x) format(round(x, 1), nsmall = 1)) + 
  scale_color_gradient("Scaled incidence", low = "white", high = "#134471",
                      labels = function(x) format(round(x, 1), nsmall = 1)) +
  scale_x_continuous("Time (weeks)", expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  theme(axis.title.y = element_blank(),
        legend.position = "top",
        strip.background = element_blank())

# Plot sncf
sncf_pl <- sncf_pl0 %>%
  ggplot(aes(x = x , y = y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yhig), alpha = 1, fill = "grey90") + 
  geom_line(color = "black") + 
  geom_hline(aes(yintercept = 0), color = "grey") +
  scale_x_continuous("Distance (km)", expand = c(0,0), limits = c(0,1000)) + 
  scale_y_continuous("Synchrony", limits = c(-0.1,1), expand = c(0,0)) +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank())

pl_ss3 <- sncf_sim_lag$data %>%
  ggplot(aes(y = lag, x = dis_km, color = lag)) +
  geom_ribbon(data = sncf_sim_lag$pred,
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5), 
              color = "transparent", alpha = 1, fill = "grey90") +
  geom_point(size = 3) +
  geom_smooth(data = sncf_sim_lag$pred,
              aes(y = Estimate, color = Estimate), color = "grey", linewidth = 0.7) +
  #geom_line(aes(y = lag, x = lag * (sncf_sim_lag$speed$`mean(speed)`/4))) +
  #scale_color_viridis_c(option = "inferno", direction = -1, end = 0.8, begin = 0.2) +
  scale_colour_gradient2(low = "#F2B336", mid = "white", high = "#134471", midpoint = 16) + 
  scale_y_continuous("Lag (week)") +
  scale_x_continuous("Distance (km)") +
  facet_grid(~glue(" ")) +
  theme(legend.position = "none",
        strip.background = element_blank())

pl_ss4 <- ggplot() +
  geom_polygon(data = filter(map_data("world"), region == country_name),
               aes(x=long, y = lat, group = group), fill = "grey90") +
  geom_point(data = sncf_sim_lag$data,
             aes(x = lon, y = lat, color = lag),
             size = 3) +
  #scale_color_viridis_c("Lag \n(week)", option = "inferno", direction = -1, end = 0.8, begin = 0.2) +
  scale_colour_gradient2("Lag \n(week)", low = "#F2B336", mid = "white", high = "#134471", midpoint = 16) + 
  scale_x_continuous("Longitude (°)") +
  scale_y_continuous("Latitude (°N)") +
  facet_grid(~glue("{country_name}")) +
  coord_map() +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        strip.background = element_blank())

pl_sim_sncf <- (pl_ss4| sim_pl + 
                  theme(legend.justification = "center")) / (pl_ss3 | sncf_pl) 

pl_sim_sncf

ggsave(plot = pl_sim_sncf, 
       filename = "_presentationfigures/presentation_pl2.svg", 
       width = 7, 
       height = 7)



