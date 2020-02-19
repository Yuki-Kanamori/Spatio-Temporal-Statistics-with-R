## Wikle, C. K., Zammit-Mangion, A., and Cressie, N. (2019), 
## Spatio-Temporal Statistics with R, Boca Raton, FL: Chapman & Hall/CRC
## Copyright (c) 2019 Wikle, Zammit-Mangion, Cressie
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## Modified by Kota Sawada

## Before running the following line, install "needs" from CRAN, STRbook and pforeach from github.
needs::needs(tidyverse, fields, gstat, RColorBrewer, sp, spacetime, STRbook, pforeach)

## ------------------------------------------------------------------------
data("NOAA_df_1990", package = "STRbook")
Tmax <- filter(NOAA_df_1990,       # subset the data
               proc == "Tmax",     # only max temperature
               month == 7,         # July
               year == 1993)        # year of 1993
head(Tmax)

## ------------------------------------------------------------------------
pred_grid <- expand.grid(lon = seq(-100, -80, length = 20),
                         lat = seq(32, 46, length = 20),
                         day = seq(4, 29, length = 6))
head(pred_grid)

## ------------------------------------------------------------------------
Tmax_no_14 <- filter(Tmax, !(day == 14))         # remove day 14

Tmax_July_idw <- idw(formula = z ~ 1,            # dep. variable
                     locations = ~ lon + lat + day, # inputs
                     data = Tmax_no_14,          # data set
                     newdata = pred_grid,        # prediction grid
                     idp = 5)                    # inv. dist. pow.
head(Tmax_July_idw)

## ------------------------------------------------------------------------
ggplot(Tmax_July_idw) +
    geom_tile(aes(x = lon, y = lat,
                  fill = var1.pred)) +
    fill_scale(name = "degF") +    # attach color scale
    xlab("Longitude (deg)") +           # x-axis label
    ylab("Latitude (deg)") +            # y-axis label
    facet_wrap(~ day, ncol = 3) +       # facet by day
    coord_fixed(xlim = c(-100, -80),
                ylim = c(32, 46))  +    # zoom in
    theme_bw()                          # B&W theme

## ------------------------------------------------------------------------
pred_obs_dist_mat <- rdist(select(pred_grid, lon, lat, day),
                           select(Tmax_no_14, lon, lat, day))
Wt_IDW <- function(theta, dist_mat) 1 / dist_mat ^ theta
Wtilde <- Wt_IDW(theta = 5, dist_mat = pred_obs_dist_mat)
# Wtilde[k, l] = 1 / d(k-th prediction location, l-th observation location)^theta

Wtilde_rsums <- rowSums(Wtilde)
W <- Wtilde / Wtilde_rsums # weight matrix or influence matrix

## ------------------------------------------------------------------------
z_pred_IDW <- as.numeric(W %*% Tmax_no_14$z)

## ------------------------------------------------------------------------
summary(Tmax_July_idw$var1.pred - z_pred_IDW)

## ------------------------------------------------------------------------
theta <- 0.5                       # set bandwidth
Wt_Gauss <- function(theta, dist_mat) exp(-dist_mat ^ 2 / theta)
Wtilde <- Wt_Gauss(theta = 0.5, dist_mat = pred_obs_dist_mat)
Wtilde_rsums <- rowSums(Wtilde)    # normalizing factors
W <- Wtilde / Wtilde_rsums           # normalized kernel weights
z_pred2 <- W %*% Tmax_no_14$z      # predictions

Tmax_July_Gauss <- pred_grid %>%
  mutate(var1.pred = z_pred2)

ggplot(Tmax_July_Gauss) +
  geom_tile(aes(x = lon, y = lat,
                fill = var1.pred)) +
  fill_scale(name = "degF") +    # attach color scale
  xlab("Longitude (deg)") +           # x-axis label
  ylab("Latitude (deg)") +            # y-axis label
  facet_wrap(~ day, ncol = 3) +       # facet by day
  coord_fixed(xlim = c(-100, -80),
              ylim = c(32, 46))  +    # zoom in
  theme_bw()  

## ------------------------------------------------------------------------
obs_obs_dist_mat <- rdist(select(Tmax, lon, lat, day),
                          select(Tmax, lon, lat, day))

## ------------------------------------------------------------------------
LOOCV_score <- function(Wt_fun, theta, dist_mat, Z) {
  Wtilde <- Wt_fun(theta, dist_mat)
  CV <- 0
  for(i in 1:length(Z)) {
    Wtilde2 <- Wtilde[i,-i]
    W2 <- Wtilde2 / sum(Wtilde2)
    z_pred <- W2 %*% Z[-i]
    CV[i] <- (z_pred - Z[i])^2
  }
  mean(CV)
}
# LOOCV_score <- function(Wt_fun, theta, dist_mat, Z) {
#   Wtilde <- Wt_fun(theta, dist_mat)
#   pforeach(i = seq_along(along.with = Z))({
#     Wtilde2 <- Wtilde[i,-i]
#     W2 <- Wtilde2 / sum(Wtilde2)
#     z_pred <- W2 %*% Z[-i]
#     (z_pred - Z[i]) ^ 2
#   }) %>% 
#     mean()
# }
## Tested parallel calculation but it takes longer...

## ------------------------------------------------------------------------
LOOCV_score(Wt_fun = Wt_IDW,
             theta = 5,
             dist_mat = obs_obs_dist_mat,
             Z = Tmax$z)

LOOCV_score(Wt_fun = Wt_Gauss,
            theta = 0.5,
            dist_mat = obs_obs_dist_mat,
            Z = Tmax$z)

## ------------------------------------------------------------------------
theta_IDW <- seq(4, 6, length = 21)
theta_Gauss <- seq(0.1, 2.1, length = 21)
# CV_IDW <- CV_Gauss <- 0

# for(i in seq_along(theta_IDW)) {
#   CV_IDW[i] <- LOOCV_score(Wt_fun = Wt_IDW,
#                            theta = theta_IDW[i],
#                            dist_mat = obs_obs_dist_mat,
#                            Z = Tmax$z)
# 
#   CV_Gauss[i] <- LOOCV_score(Wt_fun = Wt_Gauss,
#                              theta = theta_Gauss[i],
#                              dist_mat = obs_obs_dist_mat,
#                              Z = Tmax$z)
# }

CV_IDW <- pforeach(i = seq_along(theta_IDW))({
  LOOCV_score(Wt_fun = Wt_IDW,
              theta = theta_IDW[i],
              dist_mat = obs_obs_dist_mat,
              Z = Tmax$z)
})
CV_Gauss <- pforeach(i = seq_along(theta_Gauss))({
  LOOCV_score(Wt_fun = Wt_Gauss,
              theta = theta_Gauss[i],
              dist_mat = obs_obs_dist_mat,
              Z = Tmax$z)
})

## ------------------------------------------------------------------------
# par(mfrow = c(1,2))
# plot(theta_IDW, CV_IDW,
#      xlab = expression(alpha),
#      ylab = expression(CV[(m)](alpha)),
#      ylim = c(7.4, 8.5), type = 'o')
# plot(theta_Gauss, CV_Gauss,
#      xlab = expression(theta),
#      ylab = expression(CV[(m)](theta)),
#      ylim = c(7.4, 8.5), type = 'o')

CV <- data.frame(theta = theta_IDW, CV = CV_IDW, method = "IDW", stringsAsFactors = FALSE) %>% 
  bind_rows(data.frame(theta = theta_Gauss, CV = CV_Gauss, method = "Gauss", stringsAsFactors = FALSE))
ggplot(CV, aes(x = theta, y = CV)) +
  geom_point() +
  geom_line() +
  facet_wrap(~method, scales = "free_x") +
  theme_bw() +
  labs(x = expression(paste("Bandwidth ", theta, " / ", alpha)), y = expression(CV[(m)]))

## ------------------------------------------------------------------------
theta_IDW[which.min(CV_IDW)]
min(CV_IDW)

## ------------------------------------------------------------------------
theta_Gauss[which.min(CV_Gauss)]
min(CV_Gauss)
