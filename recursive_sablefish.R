library(nwfscSurvey)
library(sdmTMB)
library(tidyverse)
library(ggplot2)

bio <- pull_bio(common_name="sablefish", survey="NWFSC.Combo")
names(bio) <- tolower(names(bio))
haul <- pull_haul(survey="NWFSC.Combo")
names(haul) <- tolower(names(haul))

# filter ages to be in core range
d <- bio
d <- dplyr::filter(d, sex == "F", !is.na(age))

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
  dplyr::filter(age <= 6)

subset <- dplyr::left_join(dplyr::select(subset, -lat, -lon,-year), locs) |>
  dplyr::left_join(yrs)
subset <- add_utm_columns(subset, ll_names = c("lon", "lat"))
subset$fyear <- as.factor(subset$year)
subset$fage <- as.factor(subset$age)
subset$cohort <- (subset$year - subset$age)
subset$present[which(is.na(subset$present))] <- 0

max_year <- max(subset$year)
glms <- list()

for(a in min(subset$age):(max(subset$age) - 1)) {
  # Take all age 'a' fish from 2003 - 2022. construct a coarse mesh (small n)
  subset_age <- dplyr::filter(subset, age == a, year < max_year, n_total > 0)
  mesh <- make_mesh(subset_age, xy_cols = c("X","Y"), cutoff=50)
  subset_age$notn <- subset_age$n_total - subset_age$n
  # Fit basic model -- include RW in intercept and spatiotemporal effects (missing 2020)
  fit <- sdmTMB(cbind(n, notn) ~ 1,
                #time_varying = ~ 1, # time-varying intercept
                           spatial = "on", # spatial field on
                           spatiotemporal="rw", # random walk in spatiotemporal
                           time="year",
                           mesh=mesh,
                           family = binomial(),
                           data=subset_age,
                extra_time = (2003:2023)[which(2003:2023 %in% as.numeric(names(table(subset_age$year))) ==FALSE)])
  
  # predict to locations in next age / year
  pred_df <- dplyr::filter(subset, age == (a+1), year > min(subset$year), n_total > 0)

  pred <- predict(fit, pred_df)
  pred$expected_n <- plogis(pred$est) * pred$n_total # predicted number of fish the next year
  glms[[a + 1]] <- glm(n ~ log(expected_n), data = pred)
}

