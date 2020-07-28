library(tidyverse)
library(mgcv)
library(sf)
library(rgdal)
library(raster)
library(lwgeom)
library(forecast)
library(magrittr)
library(tsibble)
library(geosphere)
source(here::here("R/sfc_as_cols.R"))
load(here::here("data/SouthBaySeagrass.rdata"))
load(here::here("data/HogIslandBaySeagrass.rdata"))
load(here::here("data/sample_locs.rdata"))

sg <- bind_rows(hi_sg, sb_sg)

sites <- sample_locs %>% 
  st_transform(.,st_crs(sg))


sg_dist <- NULL
for (i in 2012:2018){
  for (j in 1:length(unique(sites$site))){
    for (k in c("fall","summer")){
      message(paste(i, unique(sites$site)[j], k))
      
      site_dist <- 
        sites %>% 
        dplyr::filter(year == i, 
                      site == unique(sites$site)[j],
                      season == k) %>% 
        ungroup()
      
     inside <- sg %>%
        filter(YEAR == i) %>%
        st_union() %>%
        st_as_sf() %>%
        st_intersection(.,site_dist) %>% 
        nrow()
     
     sg_multipoint <- sg %>% 
       filter(YEAR == i) %>% 
       st_union() %>%
       st_as_sf() %>% 
       st_cast("MULTIPOINT") 
     
     if (nrow(sg_multipoint) == 0) next
     
     output <- site_dist %>% 
       mutate(sg_dist = as.numeric(st_distance(.,sg_multipoint)),
              is_inside = ifelse(inside == 1, TRUE, FALSE),
              sg_dist = ifelse(is_inside, sg_dist * -1, sg_dist))
      
     assign("sg_dist", rbind(
       sg_dist,
       output
     ))
    }
  }
}

sg_dist %<>% 
  mutate(sg_dist = 
    ifelse(str_detect(site, "HI") &
             year %in% c(2013, 2014, 2016),
           NA,
           ifelse(str_detect(site, "SB") &
                    year %in% c(2014, 2016),
                  NA,
                  sg_dist)))

ggplot(data = sg_dist %>% 
         filter(season == "summer",
                str_detect(site, "Bare"))) +
  geom_point(aes(x  = year, y = sg_dist, color = is_inside))

save(sg_dist, file = here::here("data/distance_to_meadow_edges.rdata"))

