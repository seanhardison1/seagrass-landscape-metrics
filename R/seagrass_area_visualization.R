library(tidyverse)
library(mgcv)
library(sf)
library(rgdal)
library(raster)
library(lwgeom)
library(forecast)
library(tsibble)
load(here::here("data/SouthBaySeagrass.rdata"))
load(here::here("data/HogIslandBaySeagrass.rdata"))

# Time series of total seagrass area for South Bay

sb_ts <- sb_sg %>% 
  dplyr::group_by(YEAR) %>% 
  summarise(do_union = F) %>% 
  mutate(sg_area =  as.numeric(units::set_units(st_area(.), "km^2"))) %>% 
  st_set_geometry(NULL) %>% 
  tidyr::complete(YEAR = full_seq(2001:2018, 1)) %>% 
  dplyr::rename(year = YEAR) %>% 
  mutate(meadow = "SB")

# Time series of total seagrass area for Hog Island Bay

hi_ts <- hi_sg %>% 
  dplyr::group_by(YEAR) %>% 
  summarise(do_union = F) %>% 
  mutate(sg_area = as.numeric(units::set_units(st_area(.), "km^2"))) %>% 
  st_set_geometry(NULL) %>% 
  tidyr::complete(YEAR = full_seq(2007:2018, 1)) %>% 
  dplyr::rename(year = YEAR) %>% 
  mutate(meadow = "HI")

sg_ts <- bind_rows(hi_ts, sb_ts) %>%
  tsibble(index = year, key = meadow) %>% 
  mutate(across(contains('sg_area'), 
                .fns = list(interp = ~na.interp(., linear = T))),
         interp = ifelse(is.na(sg_area), TRUE, FALSE))

save(sg_ts, file = here::here("data/seagrass_area.rdata"))

ggplot(sg_ts) +
  geom_point(aes(x = year, y = sg_area_interp, group= meadow), color = "purple") +
  geom_line(aes(x = year, y = sg_area, 
                group= meadow), color = "purple") +
  geom_point(aes(x = year, y = sg_area, 
                 color= meadow)) +
  geom_line(aes(x = year, y = sg_area, 
                color= meadow))



