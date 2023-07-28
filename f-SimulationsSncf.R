#######################################################################################################
# Obtain sncf for simulations from the models in different locations 
#######################################################################################################

fun_SimulationsSncf <- function(sim_id = 1) {
  # Args: 
  # sim_id: Id of the simulation (id)
  # Return: 
  # Data, results and plot sncf (list)
  
  # Prepare the data to wide 
  sncf_df <- bind_rows(simulations_l) %>%
    filter(state_var == "CC_obs") %>%
    filter(.id == sim_id) %>%
    mutate(inc = value/N) %>% 
    left_join(spat_dat[c("loc","lon","lat","loc_city_name")]) %>%
    select(loc, lon, lat, week_no, inc) %>%
    arrange(week_no) %>%
    pivot_wider(names_from = week_no,
                values_from = inc)
  
  # Estimate the nonparametric (cross-)correlation function
  sncf_result0 <- Sncf(x = sncf_df$lon, y = sncf_df$lat, z = sncf_df[,4:523], resamp=1000, latlon = TRUE)
  # sncf_result1 <- Sncf(x = sncf_df$lon, y = sncf_df$lat, z = sncf_df[,4:522], w = sncf_df[,5:523], resamp=1000, latlon = T)
  # sncf_result2 <- Sncf(x = sncf_df$lon, y = sncf_df$lat, z = sncf_df[,4:521], w = sncf_df[,6:523], resamp=1000, latlon = T)
  
  # Get labels for the plots
  all_params_labs <- all_parms %>%
    mutate(mod = case_when(eps == 0 ~ "SIS",
                           eps == 1 ~ "SIR"),
           alpha = 1/(52*alpha)) %>%
    filter(.id == sim_id)
  
  # Plot simulations
  sim_pl <- bind_rows(simulations_l) %>%
    filter(state_var == "CC_obs") %>%
    filter(.id == sim_id) %>%
    left_join(spat_dat[c("loc","lon","lat","loc_city_name")]) %>%
    mutate(loc_lab = glue("{loc_city_name} ({loc})")) %>%
    group_by(.id, loc) %>%
    mutate(inc = value/N,
           inc_minmax = (inc - min(inc))/ (max(inc) - min(inc))) %>%
    ggplot(aes(fill = inc_minmax, x = week_no/52, y = reorder(loc_lab, lat))) + 
    geom_tile() +
    scale_fill_viridis_c("Scaled incidence",
                         labels = function(x) format(round(x, 1), nsmall = 1)) + 
    scale_x_continuous("Year", expand = c(0,0), breaks = 0:10) + 
    scale_y_discrete(expand = c(0,0)) + 
    facet_grid(~glue("{all_params_labs$mod}: R0 = {all_params_labs$R0}, alpha = 1/52*{all_params_labs$alpha}")) +
    theme(axis.title.y = element_blank(),
          legend.position = "top",
          strip.background = element_blank())
  
  # Plot sncf
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
  # sncf_pl1 <- fun_sncf_pl(sncf_result = sncf_result1)
  # sncf_pl2 <- fun_sncf_pl(sncf_result = sncf_result2)
  
  # sncf_pl <- bind_rows(sncf_pl0, sncf_pl1, sncf_pl2, .id = "lag") %>%
    # mutate(lag = case_when(lag == 1 ~ "No lag",
    #                        lag == 2 ~ "Lag: 1 wk",
    #                        lag == 3 ~ "Lag: 2 wk"),
    #        lag = factor(lag, levels = c("No lag", "Lag: 1 wk", "Lag: 2 wk"))) %>%
  sncf_pl <- sncf_pl0 %>%
    ggplot(aes(x = x , y = y)) +
    geom_ribbon(aes(ymin = ylow, ymax = yhig), alpha = 0.3, fill = "#31688EFF") + 
    geom_line(color = "#31688EFF") + 
    #geom_hline(aes(yintercept = ymean), linetype = "dotted", color = "#31688EFF") +
    geom_hline(aes(yintercept = 0), color = "grey") +
    scale_x_continuous("Distance (km)", expand = c(0,0), limits = c(0,1000)) + 
    scale_y_continuous("Synchrony", limits = c(-0.1,1), expand = c(0,0)) +
    # facet_grid(~lag) + 
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_blank())
  
  #pl <- (sim_pl | sncf_pl) + plot_layout(widths = c(5,2.5))
  # pl <- (sim_pl / sncf_pl) + plot_layout(heights = c(2,1.5))
  return(list(data = sncf_df,
              sncf = sncf_result0,
              plot = list(sim_pl, sncf_pl)))
}

