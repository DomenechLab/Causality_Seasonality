#######################################################################################################
# Obtain sncf for simulations from the models in different locations 
#######################################################################################################

fun_SimulationsSncf <- function(data_sim_l = simulations_l, set_id = 1) {
  # Args: 
  # data_sim_l: List of simulations per location
  # set_id: Id of the simulation (id)
  # Return: 
  # Data, results and plot sncf (list)
  
  # Prepare the data to wide 
  sncf_df <- bind_rows(data_sim_l) %>%
    filter(state_var == "CC_obs") %>%
    filter(.id == set_id) %>%
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
    filter(.id == set_id)
  
  # Plot simulations
  sim_pl <- bind_rows(data_sim_l) %>%
    filter(state_var == "CC_obs") %>%
    filter(.id == set_id) %>%
    left_join(spat_dat[c("loc","lon","lat","loc_city_name")]) %>%
    mutate(loc_lab = glue("{loc_city_name} ({loc})")) %>%
    group_by(.id, loc) %>%
    mutate(inc = value/N,
           inc_minmax = (inc - min(inc))/ (max(inc) - min(inc))) %>%
    ggplot(aes(fill = inc_minmax, x = week_no, y = reorder(loc_lab, lat))) + 
    geom_tile() +
    scale_fill_viridis_c("Scaled incidence",
                         labels = function(x) format(round(x, 1), nsmall = 1)) + 
    scale_x_continuous("Time (weeks)", expand = c(0,0)) + 
    scale_y_discrete(expand = c(0,0)) + 
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
  
  sncf_pl <- sncf_pl0 %>%
    ggplot(aes(x = x , y = y)) +
    geom_ribbon(aes(ymin = ylow, ymax = yhig), alpha = 0.3, fill = "#31688EFF") + 
    geom_line(color = "#31688EFF") + 
    geom_hline(aes(yintercept = 0), color = "grey") +
    scale_x_continuous("Distance (km)", expand = c(0,0), limits = c(0,1000)) + 
    scale_y_continuous("Synchrony", limits = c(-0.1,1), expand = c(0,0)) +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_blank())

  return(list(data = sncf_df,
              sncf = sncf_result0,
              plot = list(sim_pl, sncf_pl)))
}

