#######################################################################################################
# Obtain sncf for simulations from the models in different locations 
#######################################################################################################

fun_ClimSncf <- function(clim = "Te") {
  # Args: 
  # sim_id: Id of the simulation (id)
  # Return: 
  # Data, results and plot sncf (list)
  
  clim_lab <- ifelse(clim == "Te", "Temperature",
                     ifelse(clim == "RH_pred", "Relative Humidity",
                            ifelse(clim == "Td", "Dew temperature", "Unknown")))
  
  # Prepare the data to wide 
  sncf_df <- bind_rows(covars_l, .id = "loc") %>%
    left_join(spat_dat[c("loc","lon","lat","loc_city_name")]) %>%
    select(loc, lon, lat, week_no, all_of(clim)) %>%
    arrange(week_no) %>%
    pivot_wider(names_from = week_no,
                values_from = all_of(clim))
  
  # Estimate the nonparametric (cross-)correlation function
  sncf_result0 <- Sncf(x = sncf_df$lon, y = sncf_df$lat, z = sncf_df[,4:523], resamp=1000, latlon = T)
  # sncf_result1 <- Sncf(x = sncf_df$lon, y = sncf_df$lat, z = sncf_df[,4:522], w = sncf_df[,5:523], resamp=1000, latlon = T)
  # sncf_result2 <- Sncf(x = sncf_df$lon, y = sncf_df$lat, z = sncf_df[,4:521], w = sncf_df[,6:523], resamp=1000, latlon = T)
  
  # Plot climate
  sim_pl <- bind_rows(clim_dat_l, .id = "loc") %>%
    left_join(spat_dat[c("loc","lon","lat","loc_city_name")]) %>%
    mutate(loc_lab = glue("{loc_city_name} ({loc})")) %>%
    group_by(loc) %>%
    mutate(clim_minmax = (get(clim) - min(get(clim)))/ (max(get(clim)) - min(get(clim)))) %>%
    ggplot(aes(fill = clim_minmax, x = week_no/52, y = reorder(loc_lab, lat))) + 
    geom_tile() +
    scale_fill_viridis_c(glue("Scaled {clim_lab}"),
                         labels = function(x) format(round(x, 1), nsmall = 1)) + 
    scale_x_continuous("Year", expand = c(0,0), breaks = c(0:10)) + 
    scale_y_discrete(expand = c(0,0)) + 
    facet_grid(~glue("{clim_lab}")) + 
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
    geom_hline(aes(yintercept = 0)) +
    scale_x_continuous("Distance (km)", expand = c(0,0), limits = c(0,1000)) + 
    scale_y_continuous("Synchrony", limits = c(-0.1,1), expand = c(0,0)) +
    # facet_grid(~lag) + 
    theme(strip.background = element_blank())
  
  return(list(data = sncf_df,
              sncf = sncf_result0,
              plot = list(sim_pl, sncf_pl)))
}

