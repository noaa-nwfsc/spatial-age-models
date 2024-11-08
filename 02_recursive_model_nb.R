library(nwfscSurvey)
library(sdmTMB)
library(tidyverse)
library(ggplot2)
library(DHARMa)
library(surveyjoin)

spp_name <- c("Pacific hake", "sablefish")[1]

if(spp_name == "Pacific hake") {
  min_age <- 1 # not many age 0s consistently sampled
  max_age <- 6
  min_year <- 2007
  max_year <- 2019
} 
if(spp_name == "sablefish") {
  min_age <- 0
  max_age <- 10
  min_year <- 2003
  max_year <- 2023
}
all_years <- min_year:max_year

bio <- pull_bio(common_name = spp_name, survey="NWFSC.Combo")
names(bio) <- tolower(names(bio))
haul <- pull_haul(survey="NWFSC.Combo")
names(haul) <- tolower(names(haul))

# filter ages to be in core range
d <- bio
d <- dplyr::filter(d, sex == "F", !is.na(age), 
                   year >= min_year,
                   year <= max_year)

# group data by age and year 
grouped <- dplyr::group_by(d, age, year, trawl_id) |>
  dplyr::summarise(n = n(), # n is the number of observations of this age in this year
                   present = ifelse(n > 1, 1, 0), # not used
                   lat = latitude_dd[1], 
                   lon = longitude_dd[1]) |>
  dplyr::ungroup()
grouped_total <- dplyr::group_by(d, trawl_id) |>
  dplyr::summarise(n_total = n()) |>
  dplyr::ungroup()
grouped <- dplyr::left_join(grouped, grouped_total)

# add trawl data -- after expanding across years and ages
grid <- expand.grid(trawl_id = unique(grouped$trawl_id), 
                    age = unique(grouped$age))
yrs <- dplyr::group_by(grouped, trawl_id) |>
  dplyr::summarise(year = year[1])
grid <- dplyr::left_join(grid, yrs)
locs <- dplyr::group_by(grouped, trawl_id) |>
  dplyr::summarize(lat = lat[1], lon = lon[1])

d <- dplyr::left_join(grid[,c("trawl_id","age")], grouped, by = c("trawl_id","age")) 
d$n[which(is.na(d$n))] <- 0
d$n_total[which(is.na(d$n_total))] <- 0

# # only include trawls with aged fish , filter out only ages 0 -- 6 for modelling
subset <- dplyr::group_by(d, trawl_id) |>
  dplyr::mutate(n_total = sum(n),
    incl = ifelse(n_total > 0, 1, 0)) |>
  dplyr::filter(incl == 1) |>
  dplyr::select(-incl) |>
  dplyr::filter(age <= max_age, age >= min_age)

subset <- dplyr::left_join(dplyr::select(subset, -lat, -lon,-year), locs) |>
  dplyr::left_join(yrs)
subset <- add_utm_columns(subset, ll_names = c("lon", "lat"))
subset$fyear <- as.factor(subset$year)
subset$fage <- as.factor(subset$age)
subset$cohort <- (subset$year - subset$age)
subset$present[which(is.na(subset$present))] <- 0

max_year <- max(subset$year)
glms <- list()

survey_grid <- surveyjoin::nwfsc_grid
survey_grid <- add_utm_columns(survey_grid, ll_names = c("lon","lat"))

for(a in min(subset$age):(max(subset$age) - 1)) {
  # Take all age 'a' fish from 2003 - 2022. construct a coarse mesh (small n)
  subset_age <- dplyr::filter(subset, age == a, year < max_year, n_total > 0)
  mesh <- make_mesh(subset_age, xy_cols = c("X","Y"), cutoff=50)
  subset_age$notn <- subset_age$n_total - subset_age$n
  # Fit basic model -- include RW in intercept and spatiotemporal effects (missing 2020)
  fit <- sdmTMB(n ~ 1,
                time_varying = ~ 1, # time-varying intercept
                time_varying_type = "ar1",
                           spatial = "on", # spatial field on
                           spatiotemporal="ar1", # random walk in spatiotemporal
                           time="year",
                           mesh=mesh,
                           family = nbinom2(),
                           offset = log(subset_age$n_total),
                           data=subset_age,
                extra_time = all_years[which(all_years %in% as.numeric(names(table(subset_age$year))) ==FALSE)])
  # These plots are a simple way to make QQ plot for training data -- generally fits well
  #res <- residuals(fit, type = "mle-mvn")
  #qqnorm(res);abline(0, 1)
  
  if(a == min(subset$age)) {
    pars <- tidy(fit, "ran_pars")
    pars$age <- a
    ran_pars <- pars
  } else {
    pars <- tidy(fit, "ran_pars")
    pars$age <- a
    ran_pars <- rbind(ran_pars, pars)
  }
  # predict to locations in next age / year
  pred_df <- dplyr::filter(subset, age == (a+1), year > min(subset$year), n_total > 0)

  pred <- predict(object = fit, newdata = pred_df, offset = pred$n_total)
  pred$expected_n <- exp(pred$est) # predicted number of fish the next year
  #glms[[a + 1]] <- glm(n ~ log(expected_n), data = pred, family = "poisson")
  #sim_residuals <- simulateResiduals(fittedModel = glms[[a + 1]])
  # plot(sim_residuals) QQ plot for test data
  #pred$resid <- sim_residuals$scaledResiduals
  # ggplot(pred_df, aes(sample = resid)) +
  #   stat_qq() +
  #   stat_qq_line() +
  #   facet_wrap(~ year, scales = "free") +
  #   theme_minimal() +
  #   labs(title = "QQ Plot Faceted by Year",
  #        x = "Theoretical Quantiles", y = "Sample Quantiles")
  # 
  saveRDS(pred, paste0("predictions/",spp_name,"_",a,".rds"))
  
  
  # predict to survey grid for mapping
  #survey_grid$year <- 2019
  new_grid <- replicate_df(survey_grid, "year", all_years)
  
  pred <- predict(fit, new_grid)
  saveRDS(pred, paste0("predictions/",spp_name,"_",a,"_surveygrid.rds"))
  
}

saveRDS(ran_pars, paste0("predictions/",spp_name,"_ranpars.rds"))




hake_pars <- readRDS(paste0("predictions/","Pacific hake","_ranpars.rds"))
sablefish_pars <- readRDS(paste0("predictions/","sablefish","_ranpars.rds"))
hake_pars$species <- "Pacific hake"
sablefish_pars$species <- "sablefish"
pars <- rbind(hake_pars, sablefish_pars)

p1 <- dplyr::filter(pars, term == "range") |> 
  ggplot(aes(age, estimate, col = species)) + 
  geom_pointrange(aes(ymin=estimate-2*std.error, ymax = estimate+2*std.error), position = position_dodge(0.4)) + 
  scale_color_viridis_d(option="magma", begin=0.2, end = 0.8) +
  xlab("Age") + 
  ylab("Estimated range (km)") + 
  theme_bw() + 
  theme(legend.position = "none")

p2 <- dplyr::filter(pars, term == "sigma_O") |> 
  ggplot(aes(age, estimate, col = species)) + 
  geom_pointrange(aes(ymin=estimate-2*std.error, ymax = estimate+2*std.error), position = position_dodge(0.4)) + 
  scale_color_viridis_d(option="magma", begin=0.2, end = 0.8) +
  xlab("Age") + 
  ylab(expression(paste("Spatial ", sigma))) + 
  theme_bw() + 
  theme(legend.position = "none")

p3 <- dplyr::filter(pars, term == "sigma_E") |> 
  ggplot(aes(age, estimate, col = species)) + 
  geom_pointrange(aes(ymin=estimate-2*std.error, ymax = estimate+2*std.error), position = position_dodge(0.4)) + 
  scale_color_viridis_d(option="magma", begin=0.2, end = 0.8) +
  xlab("Age") + 
  ylab(expression(paste("Spatiotemporal ", sigma))) + 
  theme_bw() + 
  theme(legend.position = c(0.5, 0.85))

combo <- gridExtra::grid.arrange(p1, p2, p3, nrow = 1)

ggsave(combo, filename="plots/spatial_parameters.png", height = 4, width=7)
