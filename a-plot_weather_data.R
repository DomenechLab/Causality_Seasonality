#######################################################################################################
# Plot weather data in Colombia and Germany
# Check the mathematical relationship between T, Td, and RH 
#######################################################################################################

# Load packages -----------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
source("f-Pred_RH.R")
Pred_RH <- Vectorize(FUN = Pred_RH)
library(ISOweek)
theme_set(theme_bw())
par(bty = "l", las = 2)

# Load data ---------------------------------------------------------------

if(!file.exists("_data/clim_data.rds")) {
  countries_nm <- c("Colombia", "Germany", "Spain")
  dat_l <- vector(mode = "list", length = length(countries_nm))
  
  for(i in seq_along(countries_nm)) {
    files_nm <- list.files(path = sprintf("_data/_%s", countries_nm[i])) %>% 
      str_remove(pattern = ".xlsx")
    
    dat_l[[i]] <- vector(mode = "list", length = length(files_nm))
    
    for(j in seq_along(files_nm)) {
      file_nm_split <- str_split(string = files_nm[j], pattern = "_", simplify = T)
      clim_var_nm <- file_nm_split[, 1]
      loc_nm <- file_nm_split[, 2]
      
      print(sprintf("Country: %s, var: %s, location: %s", countries_nm[i], clim_var_nm, loc_nm))
      
      dat_l[[i]][[j]] <- read_xlsx(path = sprintf("_data/_%s/%s_%s.xlsx", countries_nm[i], clim_var_nm, loc_nm), 
                                   col_names = F, 
                                   col_types = c("date", "text"))
      colnames(dat_l[[i]][[j]]) <- c("day", clim_var_nm)
      dat_l[[i]][[j]][[clim_var_nm]] <- str_extract(string = dat_l[[i]][[j]][[clim_var_nm]], 
                                                    pattern = "[0-9]+.[0-9]+") %>% as.numeric()
      dat_l[[i]][[j]] <- dat_l[[i]][[j]] %>% 
        mutate(day = as.Date(day),
               loc = loc_nm, 
               country = countries_nm[i]) %>% 
        pivot_longer(cols = all_of(clim_var_nm), names_to = "clim_var", values_to = "value")
    }
    
    dat_l[[i]] <- bind_rows(dat_l[[i]])
  }
  
  dat_l <- bind_rows(dat_l)
  
  dat_l <- dat_l %>% 
    mutate(day_no = yday(day),
           week_no = isoweek(day),
           week_date = ISOweek2date(weekdate = sprintf("%s-7", ISOweek(day))),
           month_no = month(day),
           year_no = year(day)) %>% 
    select(day, day_no, week_no, week_date, month_no, year_no, country, loc, everything())
  
  # Data in wide format
  dat_wide <- pivot_wider(data = dat_l, names_from = clim_var, values_from = value) %>% 
    rename("Te" = "MeanTemperature", 
           "Td" = "MeanDewPoint", 
           "RH" = "MeanHumidity")
  
  # Some recordings have Td > T, which is not possible
  # Correct by fixing Td = T for those cases
  dat_wide <- dat_wide %>% 
    mutate(Td = if_else(Td > Te, Te, Td)) %>% 
    select(day:loc, Te, Td, RH)
  
  dat_wide <- dat_wide %>%
    mutate(RH_pred = Pred_RH(temp = Te, dewPoint = Td), 
           country = factor(country), 
           loc = factor(loc))
  
  saveRDS(object = dat_wide, file = "_data/clim_data.rds")
} else {
  dat_wide <- readRDS("_data/clim_data.rds")
}

# Plot --------------------------------------------------------------------
pl <- ggplot(data = dat_wide %>% filter(year_no == max(year_no), Te < Inf), 
             mapping = aes(x = RH, y = RH_pred, color = Te)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_viridis(option = "viridis") + 
  facet_wrap(~ country, ncol = 2) + 
  labs(color = "Te")
print(pl)

ggsave(filename = "_figures/RH.pdf", width = 8, height = 8)

vars_nm <- c("Te", "Td", "RH")

for(var_nm in vars_nm) {
  for(ct_nm in levels(dat_wide$country)) {
    
    dat_cur <- dat_wide %>% 
      filter(country == ct_nm) %>% 
      select(all_of(c("day", "loc", var_nm))) %>% 
      rename("clim_var" = var_nm)
    
    pl <- ggplot(data = dat_cur, 
                 mapping = aes(x = day, y = clim_var)) + 
      geom_line() + 
      facet_wrap(~ loc, scales = "fixed") + 
      labs(x = "Day", y = "Value", title = sprintf("Country: %s, climatic variable: %s", ct_nm, var_nm))
    print(pl)
    
    ggsave(filename = sprintf("_figures/%s-%s.pdf", ct_nm, var_nm), plot = pl, width = 12, height = 8)
  }
}

# Identify locations with similar variability in RH ---------------------------------------

# Calculate weekly averages
dat_week <- dat_wide %>% 
  filter(year_no >= 2013) %>% 
  group_by(country, loc, week_date) %>% 
  summarise(across(.cols = Te:RH, .fns = ~ mean(.x, na.rm = TRUE))) %>% 
  ungroup() %>% 
  mutate(RH_pred = Pred_RH(dewPoint = Td, temp = Te))

# Calculate coefficient of variation (CV) of every climatic variable
dat_week_CV <- dat_week %>% 
  pivot_longer(cols = Te:RH_pred, names_to = "clim_var", values_to = "value") %>% 
  group_by(country, loc, clim_var) %>%
  summarise(CV = sd(value, na.rm = T) / mean(value, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(country = factor(country), 
         loc = factor(loc), 
         clim_var = factor(clim_var))

# Plot CV across locations in Spain and Colombia
pl <- ggplot(data = dat_week_CV %>% filter(country %in% c("Spain", "Colombia", "Germany"), clim_var %in% c("Te", "RH_pred")), 
             mapping = aes(x = country, y = CV, label = loc)) + 
  geom_point() + 
  geom_text_repel(max.overlaps = Inf) + 
  facet_wrap(~ clim_var)  + 
  labs(x = "Country", y = "Coefficient of variation")
print(pl)

ggsave(filename = "_figures/CV_RH_Te.pdf", width = 12, height = 8)

# Plot
dat_week_rho <- dat_week %>% 
  group_by(country, loc) %>% 
  summarise(rho = cor(Te, RH_pred, method = "pearson", use = "complete.obs")) %>% 
  ungroup() %>% 
  mutate(country = factor(country), 
         loc = factor(loc))

pl <- ggplot(data = dat_week_rho %>% filter(country %in% c("Spain", "Colombia", "Germany")), 
               mapping = aes(x = country, y = abs(rho), label = loc)) + 
    geom_point() + 
    geom_text_repel(max.overlaps = Inf) + 
    labs(x = "Country", y = "rho(Te, RH)")
  print(pl)
ggsave(filename = "_figures/rho_RH_Te.pdf", width = 8, height = 8)


#######################################################################################################
# END
#######################################################################################################

