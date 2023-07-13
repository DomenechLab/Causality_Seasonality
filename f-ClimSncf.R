#######################################################################################################
# Obtain sncf for simulations from the models in different locations 
#######################################################################################################

fun_ClimSncf <- function(clim = "Te") {
  # Args: 
  # sim_id: Id of the simulation (id)
  # Return: 
  # Data, results and plot sncf (list)
  
  # Prepare the data to wide 
  sncf_df <- bind_rows(covars_l, .id = "loc") %>%
    left_join(spat_dat[c("loc","lon","lat","loc_city_name")]) %>%
    select(loc, lon, lat, week_no, all_of(clim)) %>%
    arrange(week_no) %>%
    pivot_wider(names_from = week_no,
                values_from = all_of(clim))
  
  # Estimate the nonparametric (cross-)correlation function
  sncf_result <- Sncf(sncf_df$lon, sncf_df$lat, sncf_df[,4:523], resamp=1000)
  
  # Plot climate
  sim_pl <- bind_rows(clim_dat_l, .id = "loc") %>%
    left_join(spat_dat[c("loc","lon","lat","loc_city_name")]) %>%
    mutate(loc_lab = glue("{loc_city_name} ({loc})")) %>%
    group_by(loc) %>%
    mutate(clim_minmax = (get(clim) - min(get(clim)))/ (max(get(clim)) - min(get(clim)))) %>%
    ggplot(aes(fill = clim_minmax, x = week_no/52, y = reorder(loc_lab, lat))) + 
    geom_tile() +
    scale_fill_viridis_c(glue("Scaled"),
                         labels = function(x) format(round(x, 1), nsmall = 1)) + 
    scale_x_continuous("Year", expand = c(0,0)) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(axis.title.y = element_blank())
  
  # Plot sncf
  sncf_pl <- data.frame(x = c(sncf_result$real$predicted$x),
                        y = c(sncf_result$real$predicted$y),
                        ylow = c(sncf_result$boot$boot.summary$predicted$y["0.025", ]),
                        yhig = c(sncf_result$boot$boot.summary$predicted$y["0.975", ])) %>%
    mutate(ymean = mean(y)) %>%
    mutate(yhig = case_when(yhig > 1 ~ 1,
                            .default = yhig),
           ylow = case_when(ylow < -0.1 ~ -0.1,
                            .default = ylow)) %>%
    ggplot(aes(x = x , y = y)) +
    geom_ribbon(aes(ymin = ylow, ymax = yhig), alpha = 0.3, fill = "#31688EFF") + 
    geom_line(color = "#31688EFF") + 
    geom_hline(aes(yintercept = ymean), linetype = "dotted", color = "#31688EFF") +
    geom_hline(aes(yintercept = 0)) +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(limits = c(-0.1,1), expand = c(0,0)) +
    ggtitle(glue("{clim}"))
  
  pl <- (sim_pl | sncf_pl) + plot_layout(widths = c(5,1))
  
  return(list(data = sncf_df,
              sncf = sncf_result, 
              plot = pl))
}

