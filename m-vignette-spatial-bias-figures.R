#######################################################################################################
# Make figures for vignette spatial bias
#######################################################################################################

# Load packages -----------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
theme_set(theme_classic() + theme(panel.grid.minor = element_blank()))
save_plot <- T # Should all the plots be saved as a pdf? 

# Load plots ---------------------------------------------------------------------------------------
# Choose between Spain and Colombia
coun_name <- "Colombia"
if(coun_name == "Spain") locref <- "LEAS" # Select reference Madrid or Gijon
if(coun_name == "Colombia") locref <- "SKRH" # Select reference Bogota or SKRH

pl_main1_inc <- readRDS(glue("_saved/_vignette_spatial_bias/vfigD_{coun_name}_simulations.rds")) + 
  theme(legend.justification = "center")
pl_main2_sncf <- readRDS(glue("_saved/_vignette_spatial_bias/vfigB_{coun_name}_sncf.rds"))
pl_main3_ccf <- readRDS(glue("_saved/_vignette_spatial_bias/vfigC_{coun_name}_ccf.rds"))
pl_main4_map <- readRDS(glue("_saved/_vignette_spatial_bias/vfigA_{coun_name}_map.rds"))
pl_main5_gp <- readRDS(glue("_saved/_vignette_spatial_bias/vfigE_{coun_name}_gpcovariance.rds")) + 
  theme(strip.background = element_blank())
pl_main6_profile <- readRDS(glue("_saved/_vignette_spatial_bias/vfigF_{coun_name}_{locref}_profile.rds")) + 
  ggtitle("Estimation with two-location transmission model") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 1, size = unit(11, "lines")))

# Plots together 
pl_up <- ((pl_main4_map| (pl_main2_sncf / pl_main3_ccf) | pl_main1_inc) + 
            plot_layout(widths = c(1,0.5,2)))

pl_vig <- pl_up / (pl_main5_gp + pl_main6_profile + plot_layout(widths = c(2,1), nrow = 1)) + 
  plot_annotation(tag_levels = "A")
print(pl_vig)

# Save figures 
if(save_plot) {
  if(coun_name == "Colombia") { 
    ggsave(pl_vig, height = 8, width = 10, 
           file = glue("_figures/Fig_04_vignette_spatial_bias_{coun_name}.pdf")) 
  }
  if(coun_name == "Spain") { 
    ggsave(pl_vig, height = 8, width = 10, 
           file = glue("_figures/Suppfig_04_vignette_spatial_bias_{coun_name}.pdf")) 
  }
}

#######################################################################################################
# End
#######################################################################################################