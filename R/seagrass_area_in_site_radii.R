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


sg_coverage <- NULL
for (d in c(500, 1000,2000,3000)){
  
  message(d)
  Sys.sleep(0.25)
  buffed_locs <- 
    sample_locs %>% 
    st_transform(.,st_crs(sg)) %>% 
    st_as_sf() %>% 
    st_buffer(.,dist = d) %>% 
    st_as_sf() %>% 
    ungroup() %>% 
    mutate(circle_area = as.numeric(st_area(.)))
  
  for (i in 2012:2018){
    for (j in 1:length(unique(buffed_locs$site))){
      for (k in c("fall","summer")){
        message(paste(i, unique(buffed_locs$site)[j], k))
        
        int <- sg %>% dplyr::filter(YEAR == i)
        
        buffed_int <- 
          buffed_locs %>% 
          dplyr::filter(year == i, 
                        site == unique(buffed_locs$site)[j],
                        season == k) %>% 
          st_intersection(.,int) %>% 
          group_by(year, season, site, circle_area) %>% 
          summarise() %>% 
          mutate(sg_area = as.numeric(st_area(.)),
                 sg_cover = sg_area/circle_area,
                 radius = d) %>% 
          ungroup() 
        
        assign("sg_coverage", bind_rows(buffed_int, sg_coverage)) 
      }
    }
  }  
}

seagrass_coverage <- 
  sg_coverage %>% 
    st_set_geometry(NULL) %>% 
    mutate(year = as.numeric(as.character(year)),
           site = as.factor(site)) %>% 
    group_by(radius) %>% 
    dplyr::select(-circle_area, -sg_area) %>% 
    pivot_wider(names_from = radius,
                names_prefix = "sg_cover_",
                values_from = sg_cover,
                names_sep = "_", values_fill = 0) 


save(seagrass_coverage, file = here::here("data/seagrass_cover.rdata"))
