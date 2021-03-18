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

syn_sg2 <- readxl::read_xlsx(here::here("data/2018 Synoptic Seagrass Density.xlsx")) %>% 
  dplyr::select(PLOT, SHOOTS = `#shoots/m2`, SiteName = PLOT) %>% 
  mutate(SIZE = NA, DENSITY = NA, SAMPYR = 2018, AGE = NA)

syn_sg %<>% 
  bind_rows(syn_sg2) %>% 
  as_tibble()

# read in DEM
dem_raw <- raster(here::here("data/V2_topobathy_tif.tif"))

# Load VCR polygon to trim the DEM and SAV data
sb <- st_read(here::here("data/SouthBay.kml")) %>% 
  sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")) %>% 
  st_zm() 

hi <- st_read(here::here("data/HogIslandBay.kml")) %>% 
  sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")) %>% 
  st_zm()

dem_poly <- bind_rows(sb, hi) %>% 
  as(., "Spatial")
vcr_poly <- dem_poly %>% as("sf")

# crop and mask dem according to VCR polygon
dem <- mask( crop( dem_raw,y = extent(dem_poly) ), dem_poly)
dem_resample <- aggregate(dem, fact = 5, fun = mean)
dem_resample[dem_resample > 0 ] <- NA

dem_df <- dem_resample %>% 
  as("SpatialPixelsDataFrame") %>%
  as_tibble() %>% 
  dplyr::rename(depth = 1,
                longitude = 2,
                latitude = 3) %>% 
  mutate(year = 2018)

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
  mutate(depth = extract(dem_resample, as_Spatial(.)))

# intersect VIMS and synoptic data
sg_bound <- NULL
for (i in c(2015, 2017, 2018)){
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
                year = SAMPYR,
                longitude = Easting__X,
                latitude = Northing__) %>% 
  mutate(shoots_raw = shoots) %>% 
  mutate_at(vars(depth, density, shoots),
            center.and.rescale) %>% 
  mutate(year = as.factor(year),
         sg_pres = factor(ifelse(vims_dens > 0 | shoots_raw > 0, 1, 0)))

save(sg_bound, file = here::here("data/bound_sg.rdata"))

ggplot(sg_bound) +
  geom_point(aes(x = shoots, y = depth)) +
  geom_smooth(aes(x = shoots, y = depth),
              method = "lm") +
  facet_wrap(.~year)


r <- raster(extent(dem_resample),
            res = res(dem_resample),
            crs = crs(dem_resample))
vims_sg_rast <- 
  vims_sg_proc %>% 
  filter(year == 2018) %>% 
  st_cast("MULTIPOLYGON") %>% 
  fasterize::fasterize(.,r, field = "vims_dens")

vims_sg_rast[is.na(dem_resample)] <- NA
plot(vims_sg_rast)

dem_df$vims_dens <- raster::extract(vims_sg_rast, y = dem_df[,c("longitude","latitude")])

save(dem_resample, dem_df, vims_sg_proc,
     vims_sg_rast,file = here::here("data/processed_dem.rdata"))
