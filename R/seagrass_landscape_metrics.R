library(tidyverse)
library(rvest)
library(sf)
library(rgdal)
library(raster)
library(landscapetools)
library(landscapemetrics)
library(magrittr)
source(here::here("R/sfc_as_cols.R"))
load(here::here("data/sample_locs.rdata"))

sg_dens <- NULL
for (i in c(2012, 2013, 2014, 2015,2017, 2018)){
    assign("sg_dens", bind_rows(sg_dens, st_read(here::here("data",i)) %>% 
                             mutate(YEAR = i) %>% 
                             filter(DENSITY != 0))) 
}


sg_complexity <- NULL
for (d in c(500, 1000,2000,3000)){
  
  message(d)
  Sys.sleep(0.25)
  buffed_locs <- 
    sample_locs %>% 
    st_transform(.,st_crs(sg_dens)) %>% 
    st_as_sf() %>% 
    st_buffer(.,dist = d) %>% 
    st_as_sf() %>% 
    ungroup() %>% 
    mutate(circle_area = as.numeric(st_area(.)))
  
  for (i in c(2012, 2013, 2014, 2015, 2017, 2018)){
    for (j in 1:length(unique(buffed_locs$site))){
      for (k in c("fall","summer")){
        message(paste(i, unique(buffed_locs$site)[j], k))
        
        int <- sg_dens %>% dplyr::filter(YEAR == i)
        
        buffed_int <- 
          buffed_locs %>% 
          dplyr::filter(year == i, 
                        site == unique(buffed_locs$site)[j],
                        season == k) %>% 
          st_intersection(.,int) %>% 
          dplyr::select(year, season, site, DENSITY)
        
        if(nrow(buffed_int) == 0) next
        
        poly_rast <- 
          buffed_locs %>% 
          dplyr::filter(year == i, 
                        site == unique(buffed_locs$site)[j],
                        season == k) %>% 
          st_difference(.,buffed_int) %>%  
          mutate(DENSITY = 0) %>% 
          dplyr::select(year, season, site, DENSITY) %>% 
          rbind(.,buffed_int) %>% 
          dplyr::select(DENSITY) %>% 
          stars::st_rasterize() %>% 
          {. ->>  sg_int_rast} %>% 
          calculate_lsm(., level = "landscape",
                        what = c("lsm_l_ent",
                                 "lsm_l_condent",
                                 "lsm_l_joinent",
                                 "lsm_l_mutinf")) %>% 
          mutate(year = i,
                 site = unique(buffed_locs$site)[j],
                 season = k,
                 radius = d) %>% 
          dplyr::select(-layer, -level, -class, -id) %>% 
          rbind(tibble(metric = "mn_sg_dens",
                       value = mean(as.numeric(sg_int_rast$DENSITY)[!is.na(as.numeric(sg_int_rast$DENSITY))]),
                       site = unique(buffed_locs$site)[j],
                       season = k,
                       radius = d,
                       year = i))
        
        assign("sg_complexity", bind_rows(poly_rast, sg_complexity)) 
      }
    }
  }  
}

sg_complexity %<>% 
  mutate(value = ifelse(str_detect(site, "HI") & 
                                              year %in% c(2013, 2014, 2016),
                                            NA,
                                            ifelse(str_detect(site, "SB") & 
                                                     year %in% c(2014, 2016),
                                                   NA, 
                                                   value)),
                         meadow = str_extract(site, "(?:(?!_).)*")) 

save(sg_complexity, file = here::here("data/seagrass_habitat_complexity.rdata"))


ggplot(data = sg_complexity) +
  geom_point(aes(x = year, y = value, color = meadow, group = site)) +
  geom_line(aes(x = year, y = value, color = meadow, group = site)) +
  facet_wrap(radius~metric)
