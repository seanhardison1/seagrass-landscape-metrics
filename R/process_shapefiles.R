library(tidyverse)
library(rvest)
library(sf)
library(rgdal)
library(raster)
library(magrittr)
source(here::here("R/sfc_as_cols.R"))
# read in seagrass shapefiles. Special treatment for 2004 and 2006 shapefiles; which are in .e00 format
sg <- NULL
for (i in c(2001:2004, 2006:2015, 2017,2018)){
  if (i == 2004){
    assign("sg", bind_rows(sg,
                           st_read(here::here("data/beds04.e00")) %>% 
                             mutate(YEAR = i) %>% 
                             dplyr::select(YEAR) ))
  } else if (i == 2006){
    assign("sg", bind_rows(sg,
                           st_read(here::here("data/beds06.e00")) %>% 
                             mutate(YEAR = 2006) %>%  
                             dplyr::filter(TNODE_ %in% c(4667, 4717)) %>% 
                             dplyr::select(YEAR)))
  } else {
    assign("sg", bind_rows(sg, st_read(here::here("data",i)) %>% 
                             mutate(YEAR = i) %>% 
                             filter(DENSITY != 0) %>% 
             dplyr::select(YEAR)))
  }
  
}

# South Bay
sb_poly <- st_read(here::here("data/SouthBay.kml")) %>% 
  st_transform(.,st_crs(sg))

# Intersect seagrass layers with South Bay polygon
sb_sg <- 
  sg %>% 
  st_intersection(.,sb_poly)

# Fix data from 2004
oh_four_sb <-
  sb_sg %>%  
  filter(YEAR == 2004) %>% 
  st_zm() %>%
  st_cast(.,"POINT") %>% 
  sfc_as_cols(.) %>% 
  st_set_geometry(NULL) %>%
  mutate(group = as.factor(substr(rownames(.), 1, 1))) %>% 
  st_as_sf(.,coords = c("x","y")) %>% 
  group_by(group) %>% 
  summarise(n = n(), do_union = F) %>% 
  dplyr::filter(n > 2) %>% 
  st_cast(.,"POLYGON") %>% 
  mutate(YEAR = 2004,
         Name = "SouthBay") %>% 
  ungroup() %>% 
  dplyr::select(-group, -n)


# Fix data for 2006
oh_six_sb <- 
  sb_sg %>%  
  filter(YEAR == 2006) %>% 
  st_zm() %>%
  st_cast(.,"POINT") %>% 
  sfc_as_cols(.) %>% 
  st_set_geometry(NULL) %>%
  mutate(group = as.factor(substr(rownames(.), 1, 1))) %>% 
  st_as_sf(.,coords = c("x","y")) %>% 
  group_by(group) %>% 
  summarise(n = n(), do_union = F) %>% 
  st_cast(.,"POLYGON") %>% 
  mutate(YEAR = 2006,
         Name = "SouthBay") %>% 
  ungroup() %>% 
  dplyr::select(-group, -n)

# Bind fixed data 
sb_sg %<>% 
  dplyr::select(-Description) %>% 
  st_zm() %>% 
  filter(!YEAR %in% c(2004, 2006)) %>% 
  bind_rows(., oh_four_sb, oh_six_sb)

save(sb_sg, file = here::here("data/SouthBaySeagrass.rdata"))

# Hog Island

# Hog Island polygon
hi_poly <- st_read(here::here("data/HogIslandBay.kml")) %>% 
  st_transform(.,st_crs(sg))

# No imagery for HI in 2004

# Imagery exists for HI in 2006, but there is no seagrass polygon.

# Intersect HI polygon with seagrass layers
hi_sg <- 
  sg %>% 
  st_intersection(.,hi_poly) %>%   
  dplyr::select(-Description) %>% 
  st_zm() 

#Save
save(hi_sg, file = here::here("data/HogIslandBaySeagrass.rdata"))
