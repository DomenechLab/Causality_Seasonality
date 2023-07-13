#######################################################################################################
# Prepare climatic data
#######################################################################################################

CreateClimData <- function(loc_nm, n_years = 10) {
  # Args: 
  # loc_nm: name of weather station for the target location (string)
  # n_years: no of years of climatic data required (numeric)
  # Return: 
  # Table of absolute and normalized values for weekly means for Te, Td, RH, and RH_pred (data frame)
  
  # Load all climatic data
  clim_dat <- readRDS("_data/clim_data.rds")
  
  # Filter to location and given no of years 
  clim_dat <- clim_dat %>% 
    filter(loc == loc_nm, year_no >= max(year_no) - n_years + 1) %>% 
    group_by(loc, week_date) %>% 
    summarise(across(.cols = Te:RH, .fns = ~ mean(.x, na.rm = TRUE))) %>% 
    ungroup()
  
  # Merge missing weeks, if any
  week_seq_full <- seq(from = min(clim_dat$week_date), to = max(clim_dat$week_date), by = "week")
  clim_dat <- right_join(x = clim_dat, 
                         y = data.frame(week_date = week_seq_full)) %>%
    # Impute missing values for mean value:
    mutate(week = week(week_date)) %>%
    group_by(week) %>%
    mutate(Te_miss = Te,
           Td_miss = Td,
           RH_miss = RH) %>%
    mutate(loc = case_when(is.na(loc) ~ loc_nm,
                           .default = loc),
           Te = case_when(is.na(Te) ~ mean(Te, na.rm = T),
                          .default = Te),
           Td = case_when(is.na(Td) ~ mean(Td, na.rm = T),
                          .default = Td), 
           RH = case_when(is.na(RH) ~ mean(RH, na.rm = T),
                          .default = RH)) %>%
    ungroup() %>%
    # Predict RH
    mutate(RH_pred = Pred_RH(dewPoint = Td, temp = Te),
           RH_pred_miss = Pred_RH(dewPoint = Td_miss, temp = Te_miss))
  
  # Calculate normalized values
  clim_dat <- clim_dat %>% 
    mutate(week_no = as.integer((week_date - min(week_date)) / 7),
           Te_norm = Te / mean(Te, na.rm = T), 
           Td_norm = Td / mean(Td, na.rm = T), 
           RH_norm = RH / mean(RH, na.rm = T), 
           RH_pred_norm = RH_pred / mean(RH_pred, na.rm = T),
           Te_norm_miss = Te_miss / mean(Te_miss, na.rm = T), 
           Td_norm_miss = Td_miss / mean(Td_miss, na.rm = T), 
           RH_norm_miss = RH_miss / mean(RH_miss, na.rm = T), 
           RH_pred_norm_miss = RH_pred_miss / mean(RH_pred_miss, na.rm = T)) %>% 
    select(loc, week_date, week_no, everything())
  
  return(clim_dat)
}