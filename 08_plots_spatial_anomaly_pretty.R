library(tidyverse)
library(ggplot2)
library(DHARMa)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)

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

####################


states <- map_data("state")
west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
################################################################################

# adjust x axis

# plot function
dist_map2 <- function(dfile=sablefish, 
                      Z = dfile$omega_s, 
                      xlim = c(-185, -117), 
                      midp = 0, 
                      trans = 'identity',
                      axis.mod = 3,
                      max_year = 2023,
                      no_ages = 5){
  # prep maps
  library(rnaturalearth)
  library(rnaturalearthdata)
  states <- map_data("state")
  west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
  
  # adjust axis
  dfile$yrs = axis.mod*(no_ages - dfile$age - 1)
  dfile$lon2 = dfile$lon - dfile$yrs
  xmin = floor( min(dfile$lon2))
  # axis label for X
  x1 = seq(-179,-119,axis.mod )
  xages = seq(min(dfile$age), no_ages-1,1)
  x2 = rep("", length(x1)-length(xages))
  labelsx = c(x2,xages)
  xbreaks = rev(rev(x1)[1:length(xages)])
  # set midpoing for color ramp
  if(is.na(midp)==T){midp = 0.1*max(dfile$p)}
  # prep axis
  
  # begin plotting
  plotx <- 
    ggplot(data=west_coast, aes(x = long, y = lat), 
           fill = grey(0.9), color = "black") +
    geom_polygon(aes(x = long, y = lat, group=group), 
                 fill = grey(0.9), color = "black") +
    geom_point(data=dfile , 
               aes(lon2, lat, color=Z), size=0.01) +
    # might need to adjust here based on the scale of the data, esp 'trans' term.
    scale_color_gradient2(name = "Anomaly",
                          low = scales::muted("blue"), 
                          mid = "white", 
                          high = scales::muted("red"), 
                          midpoint = midp) +  # legend title 
    coord_fixed(xlim=xlim,ylim=c(32,48),ratio=1.3) +
    # might need to adjust x-axes to plot all data
    scale_x_continuous(breaks = xbreaks , labels = xages) +
    #scale_x_continuous(breaks = seq(-179, -119, axis.mod ), labels = labelsx) +
    xlab("Age-class") +
    ylab("Latitude") +
    theme_bw() + 
    theme(axis.text = element_text(size=8),
          axis.title = element_text(size=8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.7,'lines')) 
  return(plotx)
}

fig_dir = paste0(getwd(),"/plots/")

dist_map2(sablefish,no_ages = 10, xlim = c(-152, -117))
ggsave( paste0(fig_dir,"sablefish-spatial-anomaly.png"), width=4, height = 2)

dist_map2(hake,no_ages = 6, xlim = c(-140, -117))
ggsave( paste0(fig_dir,"hake_spatial_anomaly.png"), width=3, height = 2)










