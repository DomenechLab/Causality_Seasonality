####################################################################################################
# Plot spatial data in Germany, Colombia, and Spain
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
    rename(loc = station_name,
           loc_city_name = station_city,
           lon = station_lon,
           lat = station_lat) %>%
    mutate(lat_km = (dsm::latlong2km(lat = lat, lon = lon, lon0 = 0, lat0 = 0))$km.n,
           lon_km = (dsm::latlong2km(lat = lat, lon = lon, lon0 = 0, lat0 = 0))$km.e)
  saveRDS(spat_dat, glue("{dirs$data}/spat_data.rds"))
} else {
  spat_dat <- readRDS(glue("{dirs$data}/spat_data.rds"))
}

# Plot spatial data --------------------------------------------------------------------------------

# Map of Germany 
pl_germany <- ggplot() +
  geom_polygon(data = filter(map_data("world"), region == "Germany"), 
               aes(x=long, y = lat, group = group), fill = "grey90") + 
  geom_point(data = filter(spat_dat, country == "Germany"),
             aes(x = lon, y = lat, color = reorder(loc, -lat))) + 
  geom_label_repel(data = filter(spat_dat, country == "Germany"),
                   aes(x = lon, y = lat, color = reorder(loc, -lat), 
                       label = glue("{loc_city_name} ({loc})"))) + 
  scale_color_viridis_d("Station", option = "turbo") + 
  scale_x_continuous("Longitude (°)") + 
  scale_y_continuous("Latitude (°N)") + 
  coord_map() + 
  theme(legend.position = "none",
        panel.grid = element_blank())
ggsave(pl_germany, file = glue("_extra/_extrafigures/_maps/Map_Germany.pdf"), width = 5.5)

pl_spain <- ggplot() +
  geom_polygon(data = filter(map_data("world"), region == "Spain"), 
               aes(x=long, y = lat, group = group), fill = "grey90") + 
  geom_point(data = filter(spat_dat, country == "Spain"),
             aes(x = lon, y = lat, color = reorder(loc, -lat))) + 
  geom_label_repel(data = filter(spat_dat, country == "Spain"),
                   aes(x = lon, y = lat, color = reorder(loc, -lat), 
                       label = glue("{loc_city_name} ({loc})"))) + 
  scale_color_viridis_d("Station", option = "turbo") + 
  scale_x_continuous("Longitude (°)") + 
  scale_y_continuous("Latitude (°N)") + 
  coord_map() + 
  theme(legend.position = "none",
        panel.grid = element_blank())
ggsave(pl_spain, file = glue("_extra/_extrafigures/_maps/Map_Spain.pdf"), width = 5.5)

# Map of Colombia 
pl_colombia <- ggplot() +
  geom_polygon(data = filter(map_data("world"), region == "Colombia"), 
               aes(x=long, y = lat, group = group), fill = "grey90") + 
  geom_point(data = filter(spat_dat, country == "Colombia"),
             aes(x = lon, y = lat, color = reorder(loc, -lat))) + 
  geom_label_repel(data = filter(spat_dat, country == "Colombia"),
                   aes(x = lon, y = lat, color = reorder(loc, -lat), 
                       label = glue("{loc_city_name} ({loc})"))) + 
  scale_color_viridis_d("Station", option = "turbo") + 
  scale_x_continuous("Longitude (°)") + 
  scale_y_continuous("Latitude (°N)") + 
  coord_map() + 
  theme(legend.position = "none",
        panel.grid = element_blank())

ggsave(pl_colombia, file = glue("_extra/_extrafigures/_maps/Map_Colombia.pdf"), width = 5.5)

####################################################################################################
# END
####################################################################################################

