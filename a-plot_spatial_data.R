####################################################################################################
# Plot spatial data in Colombia and Spain
####################################################################################################

# Load packages ------------------------------------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
theme_set(theme_bw())
library(glue)
library(maps)
library(ggrepel)

# Set directories 
dirs <- list()
dirs$data <- "_data"
dirs$outputs <- "_outputs"
dirs$figures <- "_figures"

# Load data ----------------------------------------------------------------------------------------

if(!file.exists(glue("{dirs$data}/spat_data.rds"))) {
  spat_dat <- readxl::read_excel(glue("{dirs$data}/countries_coordinates.xlsx")) %>%
    separate(loc_city_station, sep = " ", into = c("loc", "rest")) %>%
    select(-c(rest, loc_city_airport)) %>%
    rename(loc_city_name = 3,
           lon = 4,
           lat = 5)
  
  # Load names of the locations 
  names_l <- unique(filter(readRDS("_data/clim_data.rds"), country %in% c("Colombia", "Spain"))$loc)
  names_l <- str_subset(names_l, "GCLP", negate = T)
  
  # Join the spatial data with the names of locations and select first values (main cities)
  spat_dat <- data.frame(loc = names_l) %>%
    left_join(spat_dat) %>%
    group_by(loc) %>%
    mutate(lon = mean(lon),
           lat = mean(lat)) %>%
    mutate(n = row_number()) %>%
    filter(n == 1)
  
  saveRDS(spat_dat, glue("{dirs$data}/spat_data.rds"))
} else {
  spat_dat <- readRDS("_data/spat_data.rds")
}

# Plot spatial data --------------------------------------------------------------------------------

# Map of Spain 
pl_spain <- ggplot() +
  geom_polygon(data = filter(map_data("world"), region == "Spain"), 
               aes(x=long, y = lat, group = group), fill = "grey90") + 
  geom_point(data = filter(spat_dat, loc_country_name == "Spain"),
             aes(x = lon, y = lat, color = reorder(loc, -lat))) + 
  geom_label_repel(data = filter(spat_dat, loc_country_name == "Spain"),
                   aes(x = lon, y = lat, color = reorder(loc, -lat), 
                       label = glue("{loc_city_name} ({loc})"))) + 
  scale_color_viridis_d("Station", option = "turbo") + 
  scale_x_continuous("Longitude (°)") + 
  scale_y_continuous("Latitude (°N)") + 
  coord_map() + 
  theme(legend.position = "none",
        panel.grid = element_blank())
ggsave(pl_spain, file = glue("_figures/spain_map.pdf"), width = 5.5)

# Map of Colombia 
pl_colombia <- ggplot() +
  geom_polygon(data = filter(map_data("world"), region == "Colombia"), 
               aes(x=long, y = lat, group = group), fill = "grey90") + 
  geom_point(data = filter(spat_dat, loc_country_name == "Colombia"),
             aes(x = lon, y = lat, color = reorder(loc, -lat))) + 
  geom_label_repel(data = filter(spat_dat, loc_country_name == "Colombia"),
                   aes(x = lon, y = lat, color = reorder(loc, -lat), 
                       label = glue("{loc_city_name} ({loc})"))) + 
  scale_color_viridis_d("Station", option = "turbo") + 
  scale_x_continuous("Longitude (°)") + 
  scale_y_continuous("Latitude (°N)") + 
  coord_map() + 
  theme(legend.position = "none",
        panel.grid = element_blank())
ggsave(pl_colombia, file = glue("_figures/colombia_map.pdf"), width = 5.5)

####################################################################################################
# END
####################################################################################################

