####################################################################################################
# Plot climatic data
####################################################################################################

fun_PlotClimData = function(loc = names_l[[1]],
                            country = coun_name,
                            clim_data_list = TRUE) {
  # Args: 
  # loc: name of weather station for the target location (string)
  # country: name of the country (string)
  # Return: 
  # Patchwork plot with 
  ## 1. Time series for climatic variables (including imputed values)
  ## 2. Predicted RH
  ## 3. Seasonality of the climatic variables and forcing of transmission rate 
  
  if(clim_data_list == TRUE) clim_dat_long <- clim_dat_long_l[[loc]]
  if(clim_data_list == TRUE) clim_dat <- clim_dat_l[[loc]]
  
  # Plot time series of climatic variables and imputated values 
  pl1 <- ggplot(data = clim_dat_long %>% 
                  filter(var %in% c("Te", "Td", "RH", "RH_pred", 
                                    "Te_miss", "Td_miss", "RH_miss", "RH_pred_miss")) %>%
                  mutate(miss = str_detect(var, "miss"),
                         var = str_remove(var, "_miss")),
                aes(x = week_date, y = value, color = miss)) + 
    geom_line() + 
    scale_color_manual(values = c("red","black")) + 
    facet_wrap(~ var, scales = "free_y", ncol = 2, strip.position = "right") + 
    theme(legend.position = "none",
          strip.background = element_blank())
  # Plot measured and predicted RH
  pl2 <- ggplot(data = clim_dat, aes(x = RH, y = RH_pred)) + 
    geom_point() + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed")
  # Plot time series of measured and predicted RH
  pl3 <- ggplot(data = clim_dat, aes(x = week_date, y = RH)) + 
    geom_line() + 
    geom_line(aes(x = week_date, y = RH_pred), color = "red")
  # Season plots of all climatic variables
  pl4 <- ggplot(data = clim_dat_long %>% 
                  filter(var %in% c("Te_norm", "Td_norm", "RH_pred_norm")), 
                aes(x = isoweek(week_date), y = value, 
                    group = factor(year(week_date)))) + 
    geom_line(color = "grey") + 
    geom_smooth(aes(x = isoweek(week_date), y = value, group = NULL, color = NULL)) + 
    facet_wrap(var ~ ., scales = "fixed", ncol = 3, strip.position = "top") + 
    theme(strip.background = element_blank())
  # Plot seasonal forcing of transmission rate
  pl5 <- ggplot(data = clim_dat, aes(x = week_date, y = beta_seas)) + 
    geom_line()
  pl_climate <- pl1 / (pl2 | pl3) / (pl4 | pl5) + plot_annotation(title = loc)
  ggsave(pl_climate, file = glue("_figures/_climate/Climate_{country}_{loc}.pdf"), height = 5.5)
  
  return(pl_climate)
}
