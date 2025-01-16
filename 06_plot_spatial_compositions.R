library(tidyverse)
library(ggplot2)

spp_name <- c("Pacific hake", "sablefish")[1]

if(spp_name == "Pacific hake") {
  min_age <- 1 # not many age 0s consistently sampled
  max_age <- 6
  years <- 2010:2019
  ages <- seq(min_age, min_age+4)
} 
if(spp_name == "sablefish") {
  min_age <- 0
  max_age <- 10
  years <- 2014:2023
  ages <- seq(min_age, min_age+4)
}

for(a in min_age:(max_age - 1)) {
  pred <- readRDS(paste0("predictions/",spp_name,"_",a,"_surveygrid.rds"))
  pred <- dplyr::select(pred, lon, lat, X, Y, year, est) |>
    dplyr::mutate(age = a) |>
    dplyr::filter(year %in% years)
  if(a == min_age) {
    pred_all <- pred
  } else {
    pred_all <- rbind(pred_all, pred)
  }
}

# now convert the 'est' to normalized probabilities across ages
pred_all <- dplyr::mutate(pred_all, p = exp(est)) |>
  #dplyr::group_by(X,Y,year) |>
  #dplyr::mutate(p_norm = p / sum(p)) |>
  dplyr::filter(age %in% ages)

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

# Plot coast with custom land and ocean colors
g <- ggplot(coast_proj) + 
  geom_tile(data = dplyr::filter(pred_all, 
                          p <= quantile(pred_all$p,0.99)), 
                          aes(x = X * 1000, y = Y * 1000, col = p)) + 
  geom_sf(fill = land_color) +  # Set land color
  #scale_color_viridis(option="magma", begin = 0.2, end = 0.8, name = "CPUE", trans="sqrt") +  # legend title 
  scale_color_gradient2(name = "CPUE", trans = "sqrt", low = "white", high = scales::muted("red")) + 
  theme_light() + 
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = ocean_color),  
        strip.text = element_text(color = "black"),  # Set facet label text color to black
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # Set facet label background to white
  facet_grid(age ~ year)

g

ggsave(g, filename = paste0("plots/spatial_composition_",spp_name,".png"), height = 7, width = 7)

# Make subsets for the paper based on a few dominant cohorts
# For hake pick fish age-1 in 2011 and 2015
# For sablefish, age-0 in 2013 / 2016 / 2020

if(spp_name == "Pacific hake") {
  pred_subset <- dplyr::mutate(pred_all, cohort = year - age) |>
    dplyr::filter(cohort %in% c(2010, 2014))
} else {
  pred_subset <- dplyr::mutate(pred_all, cohort = year - age) |>
    dplyr::filter(cohort %in% c(2013, 2016, 2019))
}


g <- ggplot(coast_proj) + 
  geom_point(data = dplyr::filter(pred_subset, p <= quantile(pred_subset$p,0.99)), aes(x = X * 1000, y = Y * 1000, col = p),size=0.05) + 
  geom_sf(fill = land_color) +  # Set land color
  scale_color_gradient2(name = "CPUE", high = scales::muted("red")) +  # legend title 
  theme_light() + 
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = ocean_color),  
        strip.text = element_text(color = "black"),  # Set facet label text color to black
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # Set facet label background to white
  facet_grid(cohort ~ age)

ggsave(g, filename = paste0("plots/spatial_composition_subset_",spp_name,".png"), height = 7, width = 7)

