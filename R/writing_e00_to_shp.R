library(tidyverse)
library(rvest)
library(sf)
library(rgdal)
library(raster)
library(magrittr)
library(plotly)
source(here::here("R/sfc_as_cols.R"))
ogrListLayers(here::here("data/beds04.e00"))

oh_four <- st_read(here::here("data/beds04.e00"),
                   layer = "LAB")

oh_four_hidden <- st_read(here::here("data/beds04.e00"),
                   layer = "ARC")

vcr_poly <- st_read(here::here("data/VCR_poly3.kml")) %>% 
  st_transform(.,st_crs(oh_four))

dens_extract <- 
  oh_four %>% 
  st_intersection(., vcr_poly) %>% 
  st_zm()

fixed <- oh_four_hidden %>% 
  st_intersection(., vcr_poly) %>% 
  st_zm() %>% 
  st_cast(.,"POINT") %>% 
  sfc_as_cols(.) %>% 
  st_set_geometry(NULL) %>%
  mutate(group = as.factor(substr(rownames(.), 1, 4))) %>% 
  st_as_sf(.,coords = c("x","y"),
           crs = crs(oh_four)) %>% 
  group_by(group) %>% 
  summarise(n = n(), do_union = F) %>% 
  dplyr::filter(n > 4) %>% 
  st_cast(.,"POLYGON") %>% 
  st_join(.,dens_extract) %>% 
  filter(!is.na(DENSITY)) %>% 
  ungroup() %>% 
  dplyr::select(-group, -n)

st_write(fixed, here::here("data/beds04.shp"), driver="ESRI Shapefile",
         append = F)  

