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
load(here::here("data/sample_locs.rdata"))


sg <- bind_rows(hi_sg, sb_sg)

buffed_locs <- 
  sample_locs %>% 
  st_transform(.,st_crs(sg)) %>% 
  st_as_sf() %>% 
  st_buffer(.,dist = 1000) %>% 
  st_as_sf() %>% 
  ungroup()

t <- buffed_locs %>% 
  filter(year == 2012, season == "summer") %>% 
  dplyr::select(insideSeagrass)


ggplot() +
  geom_sf(data = t)


out <- NULL
for (i in 2012:2018){
  int <- sg %>% dplyr::filter(YEAR == i)
  
  buffed_int <- 
    buffed_locs %>% 
    dplyr::filter(year == i, season == "summer") %>% 
    st_intersection(.,int) %>% 
    group_by(site) %>% 
    summarise(do_union = F) %>% 
    mutate(sg_area = st_area(., by_feature = F),
           year = i) %>% 
    ungroup()

    assign("out", bind_rows(buffed_int, out))
}  


library(plotly)
ggplot() +
  geom_sf(data = int) +
  geom_sf(data = t, color = "red") +
  geom_sf(data = buffed_int, fill = "transparent", color = "blue") 
  

plot(buffed_int %>% dplyr::select(site))
plot(int, add = T)

ggplot() +
  geom_sf(data = out) +
  facet_wrap(.~year)