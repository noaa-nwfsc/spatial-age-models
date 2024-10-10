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
  
  output <- data.frame(age = seq(min_age, max_age), mean = NA, lower = NA, upper = NA)
  
  for(a in min_age:(max_age - 1)) {
    d <- readRDS(paste0("predictions/",spp_name,"_",a,".rds"))
    g <- glm(n ~ log(expected_n), data = d, family = "poisson")
    output[which(output$age==a),"mean"] <- coef(g)[2]
    output[which(output$age==a),c("lower","upper")] <- confint(g)[2,]
  }
  
  output$species <- spp_name
  output <- dplyr::filter(output, !is.na(mean))
  if(spp == 1) {
    all_output <- output
  } else {
    all_output <- rbind(all_output, output)
  }
}

ggplot(all_output, aes(age, mean, col = species)) + 
  geom_hline(aes(yintercept=0), col="red", alpha=0.3) + 
  geom_pointrange(aes(ymin = lower, ymax = upper), position=position_dodge(0.4)) + 
  geom_point(position=position_dodge(0.4)) + 
  scale_color_viridis_d(option="magma", begin=0.2, end = 0.8) +
  xlab("Age") + 
  ylab("Estimated slope") + 
  coord_flip() + theme_bw()

ggsave("plots/glm_coefficients.png", width=7)



