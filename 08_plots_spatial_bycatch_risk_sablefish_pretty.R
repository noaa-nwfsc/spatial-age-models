library(tidyverse)
library(ggplot2)
library(DHARMa)
library(viridis)
home_dir = getwd()

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

sablefish <- dplyr::filter(d,year %in% seq(2015,2023,by=2)) |>
  dplyr::mutate(p = exp(est))

axis.mod = 3
max_year = 2023
sablefish$yrs = axis.mod*(sablefish$year-max_year)
sablefish$lon2 = sablefish$lon + sablefish$yrs
xmin = floor( min(sablefish$lon2))

# background map info ##########################################################
states <- map_data("state")
west_coast <- subset(states, region %in% c("california", "oregon", "washington","idaho","nevada","montana"))
################################################################################

# plot function
dist_maps <- function(dfile, 
                      xlim = c(-185, -117), 
                      midp = NA, 
                      trans = 'identity',
                      scalebiomass = T) {
  # data preparation
  # dfile = pred_all[pred_all$age==1,]
  dfile = dplyr::filter(dfile, p <= quantile(dfile$p,0.99))
  
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
                 fill = "wheat", color = "grey70") +
    geom_point(data=dfile , 
               aes(lon2, lat, color=p), size=0.01) +
    # might need to adjust here based on the scale of the data, esp 'trans' term.
    # scale_color_viridis_c(option = 'turbo', direction = 1,  
    #                     name='CPUE') +
    scale_color_gradient2(name = "CPUE age-1", 
                          transform = trans, 
                          low = scales::muted("black"),
                          midpoint = midp,
                          #high = 'red')+
                          high = scales::muted("red")) + 
    coord_fixed(xlim=xlim,ylim=c(32,48),ratio=1.3) +
    # might need to adjust x-axes to plot all data
    scale_x_continuous(breaks = seq(-179,-119, axis.mod ), 
                       minor_breaks = -117:-179, labels = 2003:max_year) +
    #scale_y_continuous(breaks = c(35,40,45), minor_breaks = c(32:48)) +
    xlab("Year") +
    ylab("Latitude") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.size = unit(1,'lines')) 
  return(plotx)
}

p1 = dist_maps(sablefish, xlim=c(-150,-117))



# read in ports
ports_rad <- readRDS( paste0(home_dir, "/ports_rad.rds"))

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

cutoff <- 282
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
  scale_x_continuous(breaks= seq(2003,2023,5), minor_breaks = 2003:2023) + 
  scale_color_viridis_d(option="magma", begin=0.2, end=0.8, name = "Port")


graphics.off()
png("plots/sablefish_spatial_risk.png", width=5, heigh=5, res=300, units = 'in')
patchwork::plot_layout(p1 / p2)
dev.off()

jpeg("plots/final_plots/06_sablefish_spatial_risk.jpg", width=7, height=7, res=300, units = 'in')
patchwork::plot_layout(p1 / p2)
dev.off()

