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
  mutate(sg_area = as.numeric(st_area(.))) %>% 
  st_set_geometry(NULL) %>% 
  tidyr::complete(YEAR = full_seq(2001:2018, 1)) %>% 
  dplyr::rename(year = YEAR) %>% 
  mutate(bay = "South Bay")

# Time series of total seagrass area for Hog Island Bay

hi_ts <- hi_sg %>% 
  dplyr::group_by(YEAR) %>% 
  summarise(do_union = F) %>% 
  mutate(sg_area = as.numeric(st_area(.))) %>% 
  st_set_geometry(NULL) %>% 
  tidyr::complete(YEAR = full_seq(2007:2018, 1)) %>% 
  dplyr::rename(year = YEAR) %>% 
  mutate(bay = "Hog Island Bay")

sg_ts <- bind_rows(hi_ts, sb_ts) %>% 
  tsibble(index = year, key = bay)


ggplot(sg_ts) +
  geom_point(aes(x = year, y = sg_area, color = bay)) +
  geom_line(aes(x = year, y = sg_area, color = bay))
