library(tidyverse)
library(sf)
library(raster)
library(mgcv)

# Center and rescale function
center.and.rescale <- function(x) {
  x <- x - mean(x, na.rm = T)
  x <- x / (1 * sd(x))
  return(x)
}

# read in synoptic seagrass sampling data and point locations
syn_locs <- st_read(here::here("data/Syn_Site_Points_2018")) 
syn_sg <- read.csv(here::here("data/seagrass_dens_2007-2017.csv")) %>% 
  mutate(SiteName = ifelse(str_detect(PLOT, "T"),
                           str_replace_all(PLOT, " |\\.", "-"),
                           ifelse(str_detect(PLOT, "  "),
                                  str_replace_all(PLOT, "  ", " "),
                                  PLOT))) %>% 
  dplyr::select(-PLOT)

# read in DEM
dem_raw <- raster(here::here("data/V2_topobathy_tif.tif"))

# Load VCR polygon to trim the DEM and SAV data
dem_poly <- st_read(here::here("data/VCR_poly3.kml")) %>% 
  sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")) %>% 
  st_zm() %>% 
  as(., "Spatial")

# crop and mask dem according to VCR polygon
dem <- mask( crop( dem_raw,y = extent(dem_poly) ), dem_poly)

# read in SAV shapefiles
file_vec <- paste0("sav",2011:2017)
vims_sg <- NULL
for (i in file_vec){
  assign("vims_sg", 
         bind_rows(st_read(file.path("data",i))%>% 
                     mutate(year = as.factor(str_extract(i, "\\d+"))),
                   vims_sg))
}

# intersect SAV polygons with VCR polygon
vims_sg_proc <- vims_sg %>%
  filter(DENSITY != 0) %>% 
  st_intersection(vcr_poly) %>% 
  dplyr::select(year, vims_dens = DENSITY) %>% 
  sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs"))

# join synoptic locations with shoot density data. extract depths for each sampling site
syn_sg_fin <- 
  syn_locs %>% 
  sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")) %>% 
  left_join(., syn_sg) %>% 
  filter(!is.na(DENSITY)) %>% 
  mutate(depth = extract(dem, as_Spatial(.)))

# intersect VIMS and synoptic data
sg_bound <- NULL
for (i in 2011:2017){
  int_df <- syn_sg_fin %>% 
    dplyr::filter(SAMPYR == i) %>% 
    st_join(.,vims_sg_proc %>% dplyr::filter(year == i))

  assign('sg_bound', bind_rows(int_df, sg_bound))
}

# process output
sg_bound %<>% 
  dplyr::select(-year) %>% 
  dplyr::rename(density = DENSITY,
                shoots = SHOOTS,
                year = SAMPYR) %>% 
  mutate_at(vars(depth, density, shoots),
            center.and.rescale) %>% 
  mutate(year = as.factor(year))


ggplot(sg_bound) +
  geom_point(aes(x = shoots, y = depth)) +
  geom_smooth(aes(x = shoots, y = depth),
              method = "lm") +
  facet_wrap(.~year)

# first pass at model - GAMM with random smooth by year
mod <- gam(shoots ~ s(depth) + s(year, bs = "re") + vims_dens, data = sg_bound)

