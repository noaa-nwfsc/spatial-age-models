library(tidyverse)
library(ggplot2)
library(DHARMa)
library(viridis)
library(glmmTMB)

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
    # 
    # d$age <- a
    #d$resid <- residuals(g)
    # d$resid <- qnorm(simulateResiduals(fittedModel = g)$scaledResiduals)
    # if(a == min_age) {
    #   all_d <- d
    # } else {
    #   all_d <- rbind(all_d, d)
    # }
    # Fit a GLMM with year-specific random effects
    d$x <- log(d$expected_n)
    g <- glmmTMB(n ~ x|year, data = d, family = "poisson")
    coefs <- coef(g)[[1]]$year
    coefs$species <- spp_name
    coefs$age <- a
    coefs$year <- as.numeric(rownames(coefs))
    if(a == min_age & spp == 1) {
      all_coefs <- coefs
    } else {
      all_coefs <- rbind(all_coefs, coefs)
    }
    
  }
  # all_d$species <- spp_name
  # ggplot(data = dplyr::filter(all_d_spp, species=="sablefish"), aes(sample = resid)) +
  #   stat_qq() +
  #   stat_qq_line() +  # Add a reference line
  #   labs(title = "Q-Q Plot of GLM Residuals", x = "Theoretical Quantiles", y = "Deviance Residuals") +
  #   theme_minimal()
  # 
  
  # if(spp == 1) {
  #   all_d_spp <- all_d
  # } else {
  #   all_d_spp <- rbind(all_d_spp, all_d)
  # }
  
  output$species <- spp_name
  output <- dplyr::filter(output, !is.na(mean))
  if(spp == 1) {
    all_output <- output
  } else {
    all_output <- rbind(all_output, output)
  }
}

ggplot(all_output, aes(age, exp(mean), col = species)) + 
  geom_hline(aes(yintercept=1), col="red", alpha=0.3) + 
  geom_pointrange(aes(ymin = exp(lower), ymax = exp(upper)), position=position_dodge(0.4)) + 
  geom_point(position=position_dodge(0.4)) + 
  scale_color_viridis_d(option="magma", begin=0.2, end = 0.8, name="Species") +
  xlab("Age") + 
  ylab("Estimated slope") + 
  theme_bw()

ggsave("plots/glm_coefficients.png", width=6, height = 5)


# look at coefficients from the year model, as a diagnostic
dplyr::filter(all_coefs) |>
  ggplot(aes(year, x, group=age, col = as.factor(age))) + 
  geom_line() + xlab("Year") + 
  ylab("Estimated random slope") +
  scale_color_viridis_d(option="magma", begin = 0.2, end = 0.8, name = "Age") + 
  theme_bw() + 
  theme(strip.background = element_rect(fill="white")) + 
  facet_wrap(~species, ncol = 1, scale="free_y")
ggsave("plots/glm_coefficients_time.png", width=7, height = 5)

library(tidyr)
library(ggplot2)
library(reshape2)

df_wide <- dplyr::filter(all_coefs, species == "Pacific hake") |>
  dplyr::select(x, species, age, year) |>
  pivot_wider(names_from = age, values_from = x, names_prefix = "age_")
correlation_matrix <- cor(df_wide[,-c(1:2)], use = "complete.obs") # Exclude the year column
# Melt the correlation matrix for ggplot
correlation_melt <- melt(correlation_matrix)
# Plot the heatmap
ggplot(correlation_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  labs(title = "Correlation Matrix", x = "Age", y = "Age", fill = "Correlation")



