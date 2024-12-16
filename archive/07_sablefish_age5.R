library(tidyverse)
library(ggplot2)
library(DHARMa)
library(viridis)

spp <- 2
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
    if(a == min_age) {
      all_output <- d
    } else {
      all_output <- rbind(all_output, d)
    }
  }
  
sablefish <- dplyr::filter(all_output, species == "sablefish", 
                           year %in% seq(2015,2023,by=2), age > 5) |>
  dplyr::mutate(p = exp(est)) |>
  dplyr::group_by(X, Y, year) |>
  dplyr::summarize(p = sum(p))

# sablefish <- dplyr::group_by(sablefish, year) |>
#   dplyr::mutate(centered_p = p - mean(p))

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
# g1 <- ggplot(coast_proj) + 
#   geom_point(data = sablefish, aes(x = X * 1000, y = Y * 1000, col = p), size = 0.02) + 
#   scale_color_viridis(option="magma", begin = 0.2, end = 0.8, name = "Probability") +  # legend title 
#   #scale_color_gradient2(name = "Centered \n Pr(occurrence)") + 
#   geom_sf(fill = land_color) +  # Set land color
#   theme_light() + 
#   labs(x = "Longitude", y = "Latitude") +
#   theme(panel.background = element_rect(fill = ocean_color),  
#         strip.text = element_text(color = "black"),  # Set facet label text color to black
#         strip.background = element_rect(fill = "white")) +  # Set facet label background to white
#   facet_wrap(~year, nrow = 1)

g2 <- ggplot(coast_proj) + 
  geom_point(data = sablefish, aes(x = X * 1000, y = Y * 1000, col = p), size = 0.02) + 
  scale_color_viridis(option="magma", begin = 0.2, end = 0.8, name = "CPUE") +  # legend title 
  #scale_color_gradient2(name = "Centered \n Pr(occurrence)") + 
  geom_sf(fill = land_color) +  # Set land color
  theme_light() + 
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = ocean_color),  
        strip.text = element_text(color = "black"),  # Set facet label text color to black
        strip.background = element_rect(fill = "white")) +  # Set facet label background to white
  facet_wrap(~year, nrow = 1)

p3 <- gridExtra::grid.arrange(g1, g2, nrow = 2)
ggsave(g2, filename = "plots/sablefish_spatial_age5.png", width = 7, height = 4)  

