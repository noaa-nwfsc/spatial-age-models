library(nwfscSurvey)
library(sdmTMB)

bio <- pull_bio(common_name="sablefish", survey="NWFSC.Combo")
names(bio) <- tolower(names(bio))

age_cutoff <- 10
d <- dplyr::filter(bio, age <= age_cutoff, !is.na(age), age > 0, sex!="U") %>%
  dplyr::select(trawl_id, year, datetime_utc_iso, depth_m, sex, age, latitude_dd, longitude_dd)


# create a grid of to add absences
grid <- expand.grid(age = unique(d$age), 
                    trawl_id = unique(d$trawl_id), sex = unique(d$sex))
# join in year / etc
grid <- dplyr::left_join(grid, dplyr::group_by(d, trawl_id) %>% dplyr::summarise(depth_m = depth_m[1], 
                                                                                 datetime_utc_iso = datetime_utc_iso[1], 
                                                                                 year = year[1],
                                                                                 latitude_dd = latitude_dd[1],
                                                                                 longitude_dd = longitude_dd[1]))
missing <- which(paste(grid$age, grid$trawl_id, grid$sex) %in% paste(d$age, d$trawl_id, d$sex) == FALSE)
d$present <- 1
grid$present <- 0
d <- rbind(d, grid[missing,])

# processing / setup for sdmTMB
d <- add_utm_columns(d, ll_names = c("longitude_dd", "latitude_dd"))
d$fyear <- as.factor(d$year)
d$fage <- as.factor(d$age)
d$cohort <- (d$year - d$age)
mesh <- make_mesh(d, xy_cols = c("X","Y"), cutoff=20)

# basic spatiotemporal model. 
fit <- sdmTMB(present ~ -1 + fyear,
              spatial = "on",
              spatiotemporal="iid",
              time="year",
              mesh=mesh,
              family = binomial(),
              data=d)
sanity(fit) # all good! 

# this assumes a common age spatial field, and age - specific fields

fit_age <- update(fit,
                  spatiotemporal="ar1",
                  time="age")# doesn't converge 

fit_cohort <- update(fit,
                  spatiotemporal="ar1",
                  time="cohort") # converges, 84907.66

# In these models, there's no structure by age / year / cohort -- each age-year combo
# has a unique IID spatiotemporal field. These 3 models are all doing the exact same thing:
# d$age_yr <- paste(d$age, d$year)
# fit_age_yr <- update(fit,
#                   spatiotemporal="iid",
#                   time="age_yr") # converges, AIC = 81960.48
# 
# d$cohort_yr <- paste(d$cohort, d$year)
# fit_cohort_yr <- update(fit,
#                      spatiotemporal="iid",
#                      time="cohort_yr")
# 
# d$cohort_age <- paste(d$cohort, d$age)
# fit_cohort_age <- update(fit,
#                         spatiotemporal="iid",
#                         time="cohort_age")

######### Single cohort model
d_single <- dplyr::filter(d, cohort == 2014)

mesh <- make_mesh(d_single, xy_cols = c("X","Y"), cutoff=50)

fit_single <- sdmTMB(present ~ -1 + fyear,
              spatial = "off",
              spatiotemporal="ar1",
              time="year",
              mesh=mesh,
              family = binomial(),
              data=d_single,
              extra_time = c(2020))
pred <- predict(fit_single)
ggplot(dplyr::filter(pred, year!=2020), aes(X,Y, col=est_rf)) + 
  geom_point() + 
  scale_color_gradient2() + 
  facet_wrap(~year, nrow=1)


