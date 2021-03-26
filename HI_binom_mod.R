library(tidyverse)
library(sf)
library(raster)
library(sdmTMB)
library(patchwork)
library(INLA)
library(ggtext)
library(mgcv)
library(gratia)
library(magrittr)

# Get seagrass data from presentation-------------------------
load(here::here("data/processed_sb_df.rdata"))
sg_cover_df %<>% mutate(sg_pres = ifelse(cover > 0, 1, 0)) 

# Read HI data------------------------
sb <- st_read(here::here("data/SB_INLA_poly.kml")) %>% 
  sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")) %>% 
  st_zm()

# Read and crop raster------------------------
load(here::here("data/processed_dem.rdata"))
load(here::here("data/dem_resample2.rdata"))

elevs <- 
  ggplot(dem_df) +
  geom_tile(aes(x = longitude, y = latitude, fill = elevation)) +
  scale_fill_viridis_c()

# Fit the model-----------------------

spat_mod <- gamm(sg_pres ~ elevation +
                  s(year, k = 3) +
                  te(longitude, latitude) +
                  te(longitude, latitude, year),
                 correlation = corARMA(form = ~ 1|year, p = 1),
                data = sg_cover_df,
                family = "binomial")


# evaluate
summary(spat_mod$gam)
gratia::appraise(spat_mod$gam)
draw(spat_mod$gam)

# Predict from the model-------------------
ndf <- NULL
for (i in 2011:2018){
  int <- dem_df %>% mutate(year = i)
  assign('ndf', rbind(ndf, int))
}

ndf$pred <- predict(spat_mod$gam, 
                    newdata = ndf,
                    type = "response")
save(ndf, spat_mod, file = here::here("data/sg_st_sdm_mod.rdata"))

(twed_spat_pred <- 
    ggplot(ndf) +
    geom_tile(aes(x = longitude, y = latitude, fill = pred)) +
    geom_sf(data = vims_sg_proc %>% 
              filter(vims_dens > 0), 
            fill = "transparent",
            color = "grey50") +
    scale_fill_viridis_c() +
    facet_wrap(~year) +
    labs(fill = "Probability of\n seagrass presence") +
    theme_bw() + 
    theme(axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 45),
          legend.position = c(0.85, 0.13),
          legend.title = element_text(size = 7)) +
    coord_sf(xlim = c(421000, 429000),
             ylim = c(4116000, 4127000)))


sg_df_pred <- predict(spat_mod$gam)
(twed_pred <- 
    ggplot(sg_df_pred) +
    geom_point(aes(x = elevation, y= canopy)) +
    geom_line(aes(x = elevation, y = est)) +
    theme_bw() +
    labs(y = "Shoot canopy height (cm)",
         x = "Elevation (m, NAVD88)") +
    theme(axis.title.y = element_markdown(),
          axis.title.x = element_markdown()))

load(here::here("data/bound_sg.rdata"))

sg_bound %>% 
  filter(str_detect(SiteName,"SB")) %>% 
  mutate(year = as.numeric(as.character(year))) %>% 
  group_by(year) %>% 
  dplyr::summarise(mean_shoot = mean(shoots),
                   sd_shoot = sd(shoots)) %>% 
  ggplot() +
    geom_line(aes(x = year, y = mean_shoot))
