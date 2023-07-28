#######################################################################################################
# Obtain sncf for simulations from the models in different locations 
#######################################################################################################

fun_SimulationsSncf_lag <- function(sim_id = 1, reference = "SKRH") {
  # Args: 
  # sim_id: Id of the simulation (id)
  # reference: Loc of reference (string)
  # Return: 
  # Data, results and plot sncf (list)
  
  # Prepare the data to wide
  sncf_df <- bind_rows(simulations_l) %>%
    filter(state_var == "CC_obs") %>%
    filter(.id == sim_id) %>%
    mutate(inc = value/N) %>% 
    select(loc, week_no, inc) %>%
    arrange(week_no) %>%
    pivot_wider(names_from = loc,
                values_from = inc)
  
  # Estimate the pair-wise (cross-)correlation function from reference 
  ccf_l <- list()
  for (i in names(simulations_l)) {
    ccf_l[["lag"]] <- ccf(sncf_df[[reference]], sncf_df[[i]], lag.max = 35)$lag
    ccf_l[[i]] <- ccf(sncf_df[[reference]], sncf_df[[i]], lag.max = 35)$acf
    }
  
  # Bring data together 
  ccf_max <- bind_rows(ccf_l) %>%
    pivot_longer(!lag, names_to = "loc", values_to = "ccf") %>%
    mutate(ccf = as.numeric(ccf)) %>%
    mutate(lag = as.numeric(lag)) %>%
    # Filter out negative lag
    filter(lag >= 0) %>%
    # Obtain maximum cross-correlation
    group_by(loc) %>%
    mutate(max_ccf = max(abs(ccf))) %>%
    filter(max_ccf == abs(ccf)) %>%
    ungroup() %>%
    # Obtain distance 
    left_join(spat_dat[c("loc","lon","lat","loc_city_name")]) %>%
    mutate(lat_km = (dsm::latlong2km(lat = lat, lon = lon))$km.n,
           lat_km = lat_km + abs(min(lat_km)), 
           lon_km = (dsm::latlong2km(lat = lat, lon = lon))$km.e,
           lon_km = lon_km + abs(min(lon_km))) 
  
    # Calculate distance from reference  
    lim_lat <- unlist(filter(ccf_max, loc == reference)[,"lat_km"])
    lim_lon <- unlist(filter(ccf_max, loc == reference)[,"lon_km"])

    ccf_max <- ccf_max %>%
      mutate(dis_km = sqrt((lat_km - lim_lat)^2+(lon_km  - lim_lon)^2)) 
    
    # Calculate speed 
    mod <- lm(lag ~ dis_km, data = ccf_max)
    ccf_speed <- list(mean = (1/coef(mod)[[2]])*4,
                      sd = (1/coef(mod)[[2]])*4 - (1/(coef(mod)[[2]] + sqrt(diag(vcov(mod)))[[2]]))*4) #km/mo
    
    # Plot 
    pl <- ccf_max %>%
      ggplot(aes(x = dis_km, y = lag)) +
      geom_point() + 
      geom_smooth(method = "lm") + 
      scale_x_continuous("Distance (km)") + 
      scale_y_continuous("Lag (weeks)")

  return(list(data = ccf_max,
              speed = ccf_speed,
              plot = pl))
}

