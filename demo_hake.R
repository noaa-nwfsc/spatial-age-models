library(nwfscSurvey)
library(sdmTMB)

bio <- pull_bio(common_name="Pacific hake", survey="NWFSC.Combo")
names(bio) <- tolower(names(bio))
haul <- pull_haul(survey="NWFSC.Combo")
names(haul) <- tolower(names(haul))

# truncate upper ages and filter out age 0s -- assumes age 5+ have same distribution as 5
d <- bio
d$age[which(d$age > 5)] <- 5
d <- dplyr::filter(d, age >= 1)
# Work with just females (2x samples)
d <- dplyr::filter(d, sex == "F")


# group data by age and year 
grouped <- dplyr::group_by(d, age, year, trawl_id) |>
  dplyr::summarise(n = n(), 
                   present = ifelse(n > 1, 1, 0),
                   lat = latitude_dd[1],
                   lon = longitude_dd[1])
# add trawl data -- after expanding across years and ages
grid <- expand.grid(trawl_id = unique(haul$trawl_id),
                    year = unique(grouped$year), 
                    age = unique(grouped$age))
locs <- dplyr::group_by(haul, trawl_id) |>
  dplyr::summarize(lat = latitude_dd[1], lon = longitude_dd[1])
grid <- dplyr::left_join(grid, locs)

d <- dplyr::left_join(grid, grouped) 

# processing / setup for sdmTMB
d <- add_utm_columns(d, ll_names = c("lon", "lat"))
d$fyear <- as.factor(d$year)
d$fage <- as.factor(d$age)
d$cohort <- (d$year - d$age)
d$present[which(is.na(d$present))] <- 0
mesh <- make_mesh(d, xy_cols = c("X","Y"), cutoff=50)

# basic spatiotemporal model. 
fit <- sdmTMB(present ~ -1 + fyear,
              spatial = "on",
              spatiotemporal="iid",
              time="year",
              mesh=mesh,
              family = binomial(),
              data=d)
sanity(fit) # all good! AIC 18247.06

# fit a model that allows a spatial gradient in age effects (linear)
fit_linear_svc <- sdmTMB(present ~ -1 + fyear,
              spatial = "on",
              spatiotemporal="iid",
              time="year",
              spatial_varying = ~ age,
              mesh=mesh,
              share_range = TRUE,
              family = binomial(),
              data=d)
sanity(fit_linear_svc)
AIC(fit_linear_svc) # AIC: 17178.94


# fit a model that allows different spatial effects by age
map_list = list(ln_tau_Z = factor(c(1, 1, 1, 1, 1)))
fit_svc <- sdmTMB(present ~ -1 + fyear + fage,
                         spatial = "off",
                         spatiotemporal="iid",
                         time="year",
                         spatial_varying = ~ 0 + fage,
                         mesh=mesh,
                         share_range = TRUE,
                         family = binomial(),
                         control = sdmTMBcontrol(map = map_list),
                         data=d)
sanity(fit_svc)
AIC(fit_svc) # AIC: 16926.92

# AR(1) processses in the age effects
map_list = list(ln_tau_Z = factor(rep(1, length(table(d$fyear)))))
fit_svc_age <- sdmTMB(present ~ -1 + fyear + fage,
                  spatial = "off",
                  spatiotemporal="rw",
                  time="age",
                  spatial_varying = ~ 0 + fyear,
                  mesh=mesh,
                  share_range = TRUE,
                  family = binomial(),
                  control = sdmTMBcontrol(map = map_list),
                  data=d)
sanity(fit_svc_age)
AIC(fit_svc_age) # 16894.96



# Add a smooth on the year:age fixed effect
map_list = list(ln_tau_Z = factor(rep(1, length(table(d$fyear)))))
d$year_scaled <- as.numeric(scale(d$year))
d$age_scaled <- as.numeric(scale(d$age))
fit_svc_age_smooth <- sdmTMB(present ~ 1 + s(age_scaled, year_scaled),
                      spatial = "off",
                      spatiotemporal="rw",
                      time="age",
                      spatial_varying = ~ 0 + fyear,
                      mesh=mesh,
                      share_range = TRUE,
                      family = binomial(),
                      control = sdmTMBcontrol(map = map_list),
                      data=d)
sanity(fit_svc_age_smooth)
AIC(fit_svc_age_smooth) # 16723.56

saveRDS(fit_svc_age_smooth, "sdmtmb_models/hake_model.rds")

pred <- predict(fit_svc_age_smooth)

# plot the predicted age effects
dplyr::filter(pred, year == max(d$year)) |>
  ggplot(aes(X,Y,col=epsilon_st)) + 
  facet_wrap(~ age) + 
  geom_point(size=0.05) + 
  scale_color_gradient2(midpoint=0.5) + 
  labs(col = "Age effect")
ggsave("plots/hake_spatialage_effect.png")


# plot the predicted year effects
df <- dplyr::filter(dplyr::select(pred, -year), age == max(d$age)) %>%
  pivot_longer(
    cols = starts_with("zeta_s_fyear"),  # Select the columns starting with "zeta_s_fyear"
    names_to = "year",                   # Name of the new "year" column
    names_prefix = "zeta_s_fyear",       # Remove the "zeta_s_fyear" prefix from the column names
    values_to = "zeta_s"                 # Name of the new column for the values
  )
df |>
  ggplot(aes(X,Y,col=zeta_s)) + 
  facet_wrap(~ year) + 
  geom_point(size=0.05) + 
  scale_color_gradient2() + 
  labs(col = "Year effect")
ggsave("plots/hake_spatialyear_effect.png")