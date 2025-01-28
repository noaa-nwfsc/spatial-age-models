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

spp_name <- c("Pacific hake", "sablefish")[1]

if(spp_name == "Pacific hake") {
  min_age <- 1 # not many age 0s consistently sampled
  max_age <- 6
  years <- 2003:2023
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
pred_all <- dplyr::mutate(pred_all, p = exp(est)) |>
  #dplyr::group_by(X,Y,year) |>
  #dplyr::mutate(p_norm = p / sum(p)) |>
  dplyr::filter(age %in% ages)

# NT Plot ######################################################################
axis.mod = 3
max_year = 2023
pred_all$rec_year = pred_all$year - pred_all$age
pred_all$yrs = axis.mod*(pred_all$year-max_year)
pred_all$lon2 = pred_all$lon + pred_all$yrs
xmin = floor( min(pred_all$lon2))
# pred_all$p <= quantile(pred_all$p,0.99)

# background map info ##########################################################
states <- map_data("state")
west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
################################################################################

# plot function
dist_maps <- function(dfile, 
                      xlim = c(-185, -117), 
                      midp = NA, 
                      trans = 'identity',
                      scalebiomass = T,
                      x_lab = '',
                      scale_label = 'CPUE',
                      xadjust = 0) {
  # data preparation
  # dfile = pred_all[pred_all$age==1,]
  dfile = dplyr::filter(dfile, p <= quantile(dfile$p,0.99))
  dfile$lon2 = dfile$lon2+xadjust
  # scale data to 0-1 for each year ; if decided
  
  if(scalebiomass == TRUE){
    dmean = dfile %>%
            group_by(year) %>%
            summarise(maxp = max(p, na.rm = T))
    # slow
    dfile = left_join(dfile, dmean)
    dfile$p = dfile$p/dfile$maxp
  }
  # set midpoing for color ramp
  if(is.na(midp)==T){midp = 0.1*max(dfile$p)}
  # begin plotting
  plotx <- 
    ggplot(data=west_coast, aes(x = long, y = lat), 
                 fill = grey(0.9), color = "black") +
    geom_polygon(aes(x = long, y = lat, group=group), 
               fill = grey(0.9), color = "black") +
    geom_point(data=dfile , 
               aes(lon2, lat, color=p), size=0.01) +
    # might need to adjust here based on the scale of the data, esp 'trans' term.
    # scale_color_viridis_c(option = 'turbo', direction = 1,  
    #                     name='CPUE') +
    scale_color_gradient2(name = scale_label, 
                          transform = trans, 
                          low = scales::muted("black"),
                          midpoint = midp,
                          #high = 'red')+
                          high = scales::muted("red")) + 
    coord_fixed(xlim=xlim,ylim=c(32,48),ratio=1.3) +
    # might need to adjust x-axes to plot all data
    scale_x_continuous(breaks = seq(-180,-120,axis.mod ), 
                     minor_breaks = -117:-180, labels = 2003:max_year) +
    #scale_y_continuous(breaks = c(35,40,45), minor_breaks = c(32:48)) +
    xlab(x_lab) +
    ylab("Latitude") +
    theme_bw() + 
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7),
        axis.title.x = element_text(size = 10),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(1,'lines')) 
    return(plotx)
}

# age-class each year ##########################################################

# print(m0)
# age 1 sablefish for figure 6 ####
# combine later for fig 6
m1 = dist_maps(pred_all[pred_all$age==1,],scale_label = 'CPUE age-1')
ggsave( paste0(fig_dir ,spp_name,"-age-1-year.png"), width = 6.2, height = 2)

# other age classes ####
m2 = dist_maps(pred_all[pred_all$age==2,],scale_label = 'CPUE age-2')
m3 = dist_maps(pred_all[pred_all$age==3,],scale_label = 'CPUE age-3')

if(spp_name=="sablefish"){
  m0 = dist_maps(pred_all[pred_all$age==0,],scale_label = 'CPUE age-0')
  m4 = dist_maps(pred_all[pred_all$age==4,],scale_label = 'CPUE age-4', x_lab="Year")
  ggarrange(m0,m1,m2,m3,m4, nrow = 5)}

if(spp_name=="Pacific hake"){
  m4 = dist_maps(pred_all[pred_all$age==4,],scale_label = 'CPUE age-4')
  m5 = dist_maps(pred_all[pred_all$age==5,],scale_label = 'CPUE age-5', x_lab='Year')
  ggarrange(m1,m2,m3,m4, m5, nrow = 5)}

ggsave( paste0(fig_dir , spp_name,"-age-class-year.png"), width = 6.2, height = 8)
# cohort abundance through time ################################################

# sablefish 
if(spp_name == "sablefish"){
  y2008 = dist_maps(pred_all[pred_all$rec_year==2008,],xlim = c(-170, -117), scale_label = 'CPUE 2008 recruits')
  y2010 = dist_maps(pred_all[pred_all$rec_year==2010,],xlim = c(-170, -117), scale_label = 'CPUE 2010 recruits')
  y2016 = dist_maps(pred_all[pred_all$rec_year==2016,],xlim = c(-170, -117), scale_label = 'CPUE 2016 recruits')
  y2021 = dist_maps(pred_all[pred_all$rec_year==2021,],xlim = c(-170, -117), scale_label = 'CPUE 2021 recruits')

  ggarrange(y2008, y2010, y2016, y2021,nrow = 4)

  ggsave( paste0(fig_dir ,spp_name,"-follow-age-class.png"), width = 6.2, height = 6.5)
}


# hake
if(spp_name == "Pacific hake"){
  # y2008 = dist_maps(pred_all[pred_all$rec_year==2008,],xlim = c(-170, -117))
  y2008 = dist_maps(pred_all[pred_all$rec_year==2008,],xlim = c(-170, -117), 
                    scale_label = 'CPUE 2008 recruits', xadjust = -3)
  y2010 = dist_maps(pred_all[pred_all$rec_year==2010,], xlim = c(-170, -117), 
                    scale_label = 'CPUE 2010 recruits', xadjust = -3)
  y2014 = dist_maps(pred_all[pred_all$rec_year==2014,], xlim = c(-170, -117), 
                    scale_label = 'CPUE 2014 recruits', xadjust = -3)
  y2016 = dist_maps(pred_all[pred_all$rec_year==2016,], xlim = c(-170, -117), 
                    scale_label = 'CPUE 2016 recruits', xadjust = -3)
  
  ggarrange(y2008, y2010, y2014,y2016, nrow = 4)
  
  ggsave( paste0(fig_dir ,spp_name,"-follow-age-class.png"), width = 6.2, height = 6.5)
}








