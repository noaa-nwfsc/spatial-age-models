library(sdmTMB)
# library(sdmTMBextra)
library(INLA)
library(sf)
library(tidyverse)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rasterVis)
library(readxl)
library(Matrix)
library(ggpubr)

home_dir = getwd()
fig_dir = paste0(home_dir,"/plots/")

spp_name <- c("Pacific hake", "sablefish")[2]

if(spp_name == "Pacific hake") {
  min_age <- 1 # not many age 0s consistently sampled
  max_age <- 6
  years <- 2003:2019
  ages <- seq(min_age, min_age+4)
} 
if(spp_name == "sablefish") {
  min_age <- 0
  max_age <- 10
  years <- 2003:2023
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
pred_all <- dplyr::mutate(pred_all, p = exp(est)/(1+exp(est))) |>
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
land_color  <- "wheat3"  # beige for land
ocean_color <- "grey80"  # grayish blue for ocean

# NT Plot ######################################################################
axis.mod = 3
max.year = 2023
pred_all$rec_year = pred_all$year - pred_all$age
pred_all$yrs = axis.mod*(pred_all$year-max.year)
pred_all$lon2 = pred_all$lon + pred_all$yrs
xmin = floor( min(pred_all$lon2))

# background map info ##########################################################
states <- map_data("state")
west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
################################################################################

# plot function
dist_maps <- function(dfile, xlim = c(-185, -117)) {

  prob_map <- ggplot(data=west_coast, aes(x = long, y = lat), 
                   fill = grey(0.9), color = "black") +
    geom_polygon(aes(x = long, y = lat, group=group), 
               fill = grey(0.9), color = "black") +
    geom_point(data=dfile, aes(lon2, lat, color=p), size=0.01) +
    # might need to adjust here based on the scale of the data, esp 'trans' term.
    scale_color_viridis_c(option = 'turbo', direction = 1,  
                        name='Probability') +
    coord_fixed(xlim=xlim,ylim=c(32,48),ratio=1.3) +
    # might need to adjust x-axes to plot all data
    scale_x_continuous(breaks = seq(-180,-120,axis.mod ), 
                     minor_breaks = -117:-180, labels = 2003:max_year) +
    #scale_y_continuous(breaks = c(35,40,45), minor_breaks = c(32:48)) +
    xlab("Year") +
    ylab("Latitude") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(1,'lines')) 
    return(prob_map)
}

m0 = dist_maps(pred_all[pred_all$age==0,])
m1 = dist_maps(pred_all[pred_all$age==1,])
m2 = dist_maps(pred_all[pred_all$age==2,])
m3 = dist_maps(pred_all[pred_all$age==3,])
m4 = dist_maps(pred_all[pred_all$age==4,])

ggarrange(m0,m1,m2,m3,m4,
          nrow = 5)

ggsave( paste0(fig_dir ,"Map-",spp_name,"-year.png"), width = 6.2, height = 8)


################################################################################

# sablefish 
if(spp_name == "sablefish"){
  y2008 = dist_maps(pred_all[pred_all$rec_year==2008,])
  y2010 = dist_maps(pred_all[pred_all$rec_year==2010,])
  y2016 = dist_maps(pred_all[pred_all$rec_year==2016,])
  y2021 = dist_maps(pred_all[pred_all$rec_year==2021,])

  ggarrange(y2008, y2010, y2016,y2021,nrow = 4)

  ggsave( paste0(fig_dir ,"Map-",spp_name,"-age-class.png"), width = 6.2, height = 8)
}


















