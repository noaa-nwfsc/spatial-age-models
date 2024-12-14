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
  
d <- readRDS(paste0("predictions/",spp_name,"_",a,"_surveygrid.rds"))

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

ggplot(distances_long, aes(year, value, group = Name, col = Name)) + 
  geom_line() + 
  ylab("Forecasted intensity of age-1 sablefish") + 
  xlab("Year") + 
  theme_bw() + 
  scale_color_viridis_d(option="magma", begin=0.2, end=0.8, name = "Port")

ggsave("plots/sablefish_spatial_risk.png", width = 7, height = 4)  

