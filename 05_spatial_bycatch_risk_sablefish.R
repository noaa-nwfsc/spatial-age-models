library(tidyverse)
library(ggplot2)
library(DHARMa)
library(viridis)

spp <- 2
  spp_name <- c("Pacific hake", "sablefish")[spp]
  
  if(spp_name == "Pacific hake") {
    min_age <- 1 # not many age 0s consistently sampled
    max_age <- 5
  } 
  if(spp_name == "sablefish") {
    min_age <- 0
    max_age <- 6
  }
  
  for(a in min_age:(max_age - 1)) {
    d <- readRDS(paste0("predictions/",spp_name,"_",a,"_surveygrid.rds"))
    d$species <- spp_name
    d$age <- a
    if(a == min_age) {
      all_output <- d
    } else {
      all_output <- rbind(all_output, d)
    }
  }
  
sablefish <- dplyr::filter(all_output, species == "sablefish", year %in% seq(2015,2023,by=2)) |>
  dplyr::mutate(p = plogis(est)) |>
  dplyr::group_by(year, X, Y) |>
  dplyr::mutate(p = p / sum(p)) |> # normalize by ages
  dplyr::ungroup() |>
  dplyr::filter(age <= 2) |>
  dplyr::group_by(year, X, Y) |>
  dplyr::summarise(p = sum(p)) 

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
land_color <- "#E0CDA9"  # beige for land
ocean_color <- "#D3EAF2"  # grayish blue for ocean

# Plot coast with custom land and ocean colors
ggplot(coast_proj) + 
  geom_sf(fill = land_color) +  # Set land color
  geom_point(data = sablefish, aes(x = X * 1000, y = Y * 1000, col = p)) + 
  scale_color_viridis(option="magma", begin = 0.2, end = 0.8, name = "Occurrence") +  # legend title 
  theme_light() + 
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = ocean_color),  
        strip.text = element_text(color = "black"),  # Set facet label text color to black
        strip.background = element_rect(fill = "white")) +  # Set facet label background to white
  facet_wrap(~year, nrow = 1)
ggsave("plots/sablefish_spatial_risk.png")  

