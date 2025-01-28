library(tidyverse)
library(ggplot2)
library(DHARMa)
library(viridis)

for(spp in 1:2) {
  spp_name <- c("Pacific hake", "sablefish")[spp]
  
  if(spp_name == "Pacific hake") {
    min_age <- 1 # not many age 0s consistently sampled
    max_age <- 6
  } 
  if(spp_name == "sablefish") {
    min_age <- 0
    max_age <- 10
  }
  
  for(a in min_age:(max_age - 1)) {
    d <- readRDS(paste0("predictions/",spp_name,"_",a,"_surveygrid.rds"))
    d$species <- spp_name
    d$age <- a
    if(spp == 1) {
      if(a == min_age) {
        all_output <- d
      } else {
        all_output <- rbind(all_output, d)
      }
    } else {
      all_output <- rbind(all_output, d)
    }
  }
  
}

all_output <- dplyr::filter(all_output, year == 2019)

hake <- dplyr::filter(all_output, species == "Pacific hake") |>
  dplyr::group_by(age) |>
  dplyr::mutate(omega_s = scale(omega_s))  # normalize for plotting

sablefish <- dplyr::filter(all_output, species == "sablefish") |>
  dplyr::group_by(age) |>
  dplyr::mutate(omega_s = scale(omega_s))  # normalize for plotting

map_data <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf", country = "united states of america")
# Crop the polygon for plotting and efficiency:
# st_bbox(map_data) # find the rough coordinates
coast <- suppressWarnings(suppressMessages(
  sf::st_crop(map_data,
              c(xmin = -132, ymin = 30, xmax = -117, ymax = 50))))

utm_zone10 <- 3157
coast_proj <- sf::st_transform(coast, crs = utm_zone10)

# Define the new colors
land_color <- "wheat3"  # beige for land
ocean_color <- "grey80"  # grayish blue for ocean

library(scales)
# Plot coast with custom land and ocean colors
ggplot(coast_proj) + 
  geom_point(data = hake, aes(x = X * 1000, y = Y * 1000, col = omega_s), size=0.02) + 
  geom_sf(fill = land_color) +  # Set land color
  scale_color_gradient2(name = "Anomaly",low = scales::muted("blue"), mid = "white", high = scales::muted("red"), midpoint = 0) +  # legend title 
  theme_light() + 
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = ocean_color),  
        strip.text = element_text(color = "black"),  # Set facet label text color to black
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # Set facet label background to white
  facet_wrap(~age, nrow = 1)
ggsave("plots/hake_spatial_anomalies-facet.png", width = 7, height = 5)  
  

ggplot(coast_proj) + 
  geom_point(data = sablefish, aes(x = X * 1000, y = Y * 1000, col = omega_s), size=0.02) + 
  geom_sf(fill = land_color) +  # Set land color
  scale_color_gradient2(name = "Anomaly",low = scales::muted("blue"), mid = "white", high = scales::muted("red"), midpoint = 0) +  # legend title 
  theme_light() + 
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = ocean_color),  
        strip.text = element_text(color = "black"),  # Set facet label text color to black
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # Set facet label background to white
  facet_wrap(~age, nrow = 2)
ggsave("plots/sablefish_spatial_anomalies-facet.png", width = 7, height = 8)  

