library(tidyverse)
library(sf)
library(raster)
library(magrittr)

# read in DEM
dem_raw <- raster(here::here("data/V2_topobathy_tif.tif"))

# Load VCR polygon to trim the DEM and SAV data
sb <- st_read(here::here("data/SB_INLA_poly.kml")) %>% 
  sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")) %>% 
  st_zm()

dem_poly <- sb %>% 
  as(., "Spatial")
vcr_poly <- dem_poly %>% as("sf")

# crop and mask dem according to VCR polygon
dem <- mask( crop( dem_raw,y = extent(dem_poly) ), dem_poly)
dem_resample <- aggregate(dem, fact = 50, fun = mean)
dem_resample[dem_resample > 0 ] <- NA

dem_df <- dem_resample %>% 
  as("SpatialPixelsDataFrame") %>%
  as_tibble() %>% 
  dplyr::rename(elevation = 1,
                longitude = 2,
                latitude = 3) 

save(dem_resample, dem_df, file = here::here("data/dem_resample2.rdata"))

# read in SAV shapefiles
file_vec <- paste0("sav",2011:2018)
vims_sg <- NULL
for (i in file_vec){
  assign("vims_sg", 
         bind_rows(st_read(file.path("data",i))%>% 
                     mutate(year = as.factor(str_extract(i, "\\d+"))),
                   vims_sg))
}

# intersect SAV polygons with VCR polygon
vims_sg_proc <- vims_sg %>%
  sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")) %>% 
  # filter(DENSITY != 0) %>% 
  st_intersection(vcr_poly) %>% 
  dplyr::select(year, vims_dens = DENSITY) 

# join synoptic locations with shoot density data. extract depths for each sampling site
syn_sg_fin <- 
  syn_locs %>% 
  sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")) %>% 
  left_join(., syn_sg) %>% 
  dplyr::rename(sampyr = SAMPYR) %>% 
  left_join(.,sg_bmass_canhei) %>% 
  mutate(elevation = extract(dem_resample, as_Spatial(.)))

r <- raster(extent(dem_resample),
            res = res(dem_resample),
            crs = crs(dem_resample))

sb_sg_cover <- 
  vims_sg_proc %>%
  st_cast("MULTIPOLYGON") %>% 
  group_by(year) %>% 
  fasterize::fasterize(.,r, field = "vims_dens")
  # mutate(cover = extract(dem_resample, as_Spatial(.)))

sg_cover_df <- tibble()
for (i in 2011:2018){
  if (i %in% c(2012, 2014)){
    int <- tibble(year = NA,
                  longitude = NA,
                  latitude = NA,
                  cover = NA)
  } else {
    int <- 
      vims_sg_proc %>%
      filter(year == i) %>% 
      st_cast("MULTIPOLYGON") %>% 
      fasterize::fasterize(.,r, field = "vims_dens") %>% 
      as("SpatialPixelsDataFrame") %>%
      as_tibble() %>% 
      mutate(year = i) %>% 
      dplyr::rename(cover = layer,
                    longitude = x,
                    latitude = y)
  }
  
  
  assign("sg_cover_df", rbind(sg_cover_df, int))
  
  
}


sg_cover_df$elevation <- 
  raster::extract(dem_resample, y = sg_cover_df[,c("longitude","latitude")])

save(sg_cover_df,
     file = here::here("data/processed_sb_df.rdata"))
