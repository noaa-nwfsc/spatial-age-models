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
  
d <- readRDS(paste0("predictions/",spp_name,"_",min_age+1,"_surveygrid.rds"))

sablefish <- dplyr::filter(d, 
                           year %in% seq(2015,2023,by=2)) |>
  dplyr::mutate(p = exp(est))
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

p1 <- ggplot(coast_proj) + 
  geom_point(data = dplyr::filter(sablefish, p <= quantile(sablefish$p,0.99)), aes(x = X * 1000, y = Y * 1000, col = p), size = 0.02) + 
  #scale_color_viridis(option="magma", begin = 0.2, end = 0.8, name = "CPUE", trans="sqrt") +  # legend title 
  scale_color_gradient2(name = "CPUE", low="white", high=scales::muted("red"), trans="sqrt") + 
  geom_sf(fill = land_color) +  # Set land color
  theme_light() + 
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = ocean_color),  
        strip.text = element_text(color = "black"),  # Set facet label text color to black
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +  # Set facet label background to white
  facet_wrap(~year, nrow = 1)

# read in ports
ports_rad <- readRDS("~/Documents/Github Projects/spatial-age-models/ports_rad.rds")
d$dist_AST <- sqrt((d$X - ports_rad$X[which(ports_rad$Pcid=="AST")])^2 + 
                     (d$Y - ports_rad$Y[which(ports_rad$Pcid=="AST")])^2)
d$dist_COS <- sqrt((d$X - ports_rad$X[which(ports_rad$Pcid=="COS")])^2 + 
                     (d$Y - ports_rad$Y[which(ports_rad$Pcid=="COS")])^2)
d$dist_BRK <- sqrt((d$X - ports_rad$X[which(ports_rad$Pcid=="BRK")])^2 + 
                     (d$Y - ports_rad$Y[which(ports_rad$Pcid=="BRK")])^2)
d$dist_CRS <- sqrt((d$X - ports_rad$X[which(ports_rad$Pcid=="CRS")])^2 + 
                     (d$Y - ports_rad$Y[which(ports_rad$Pcid=="CRS")])^2)
d$dist_ERK <- sqrt((d$X - ports_rad$X[which(ports_rad$Pcid=="ERK")])^2 + 
                     (d$Y - ports_rad$Y[which(ports_rad$Pcid=="ERK")])^2)
d$dist_BRG <- sqrt((d$X - ports_rad$X[which(ports_rad$Pcid=="BRG")])^2 + 
                     (d$Y - ports_rad$Y[which(ports_rad$Pcid=="BRG")])^2)
d$dist_MRO <- sqrt((d$X - ports_rad$X[which(ports_rad$Pcid=="MRO")])^2 + 
                     (d$Y - ports_rad$Y[which(ports_rad$Pcid=="MRO")])^2)

cutoff <- 232
distances <- dplyr::group_by(d, year) |>
  dplyr::summarise(AST = sum(exp(est[which(dist_AST < cutoff)])),
                   COS = sum(exp(est[which(dist_COS < cutoff)])),
                   BRK = sum(exp(est[which(dist_BRK < cutoff)])),
                   CRS = sum(exp(est[which(dist_CRS < cutoff)])),
                   ERK = sum(exp(est[which(dist_ERK < cutoff)])),
                   BRG = sum(exp(est[which(dist_BRG < cutoff)])),
                   MRO = sum(exp(est[which(dist_MRO < cutoff)])))
library(tidyr)
distances_long <- distances %>%
  pivot_longer(
    cols = -year,               # Select all columns except `year`
    names_to = "port",          # Name for the new column holding port names
    values_to = "value"         # Name for the new column holding values
  )

distances_long <- dplyr::rename(distances_long, Pcid = port)
distances_long <- dplyr::left_join(distances_long, ports_rad[,c("Pcid","Name")])
distances_long$Name[which(distances_long$Name=="Charleston (Coos Bay)")] <- "Coos Bay"

p2 <- ggplot(distances_long, aes(year, value, group = Name, col = Name)) + 
  geom_line() + 
  ylab("Forecasted CPUE of age-1 sablefish") + 
  xlab("Year") + 
  theme_bw() + 
  scale_color_viridis_d(option="magma", begin=0.2, end=0.8, name = "Port")

p3 <- gridExtra::grid.arrange(p1, p2, ncol = 1)

ggsave(p3, filename = "plots/sablefish_spatial_risk.png", width = 7, height = 6)  

