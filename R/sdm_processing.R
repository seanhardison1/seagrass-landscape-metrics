library(tidyverse)
library(sf)
library(raster)
library(magrittr)

# biomass and canopy data
load(here::here("data/sg_bmass_canhei.rdata"))

# Center and rescale function
center.and.rescale <- function(x) {
  x <- x - mean(x, na.rm = T)
  x <- x / (1 * sd(x))
  return(x)
}

# read in synoptic seagrass sampling data and point locations
syn_locs <- st_read(here::here("data/Syn_Site_Points_2018")) %>% 
  mutate(SiteName = str_replace(SiteName, "-", " "))
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
r <- raster(extent(dem),
            res = 30,
            crs = crs(dem))
# dem_resample <- aggregate(dem, fact = 10, fun = mean)
dem_resample <- resample(dem, r)
dem_resample[dem_resample > 0 ] <- NA

dem_df <- dem_resample %>% 
  as("SpatialPixelsDataFrame") %>%
  as_tibble() %>% 
  dplyr::rename(depth = 1,
                longitude = 2,
                latitude = 3) 

# read in SAV shapefiles
file_vec <- paste0("sav",2010:2018)
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

plot(vims_sg_proc %>% filter(year == 2010))
ggplot() +
  geom_sf(data = vims_sg_proc %>% filter(year == 2010)) +
  # geom_sf(data = syn_sg_fin %>% filter(sampyr == 2010), size = 0.5) +
  # geom_sf(data = int_vims) +
  geom_sf(data = int_df, aes(color = edge_dist), size = 0.5)
#Good: 2011, 2015, 2017, 2018
#Bad: 2012, 2014, 2016

#Good SB: 2013
# Bad HI: 2013

# Bad SB: 

ls_sg <- vims_sg_proc %>% 
  group_by(year) %>% 
  summarise(do.un)

# intersect VIMS and synoptic data
sg_bound <- NULL
for (i in c(2010,2011, 2013:2018)){
  
  int_vims <- vims_sg_proc %>% 
      dplyr::filter(year == i, vims_dens > 0) 
    
  if (nrow(int_vims) != 0){
      int_vims %<>% 
        summarise() %>% 
        st_cast("MULTIPOINT") 

    
      int_df <- syn_sg_fin %>% 
        dplyr::filter(sampyr == i) %>% 
        mutate(edge_dist = as.numeric(st_distance(.,int_vims))) %>% 
        st_join(.,vims_sg_proc %>% dplyr::filter(year == i)) %>% 
        mutate(edge_dist = ifelse(vims_dens == 0, edge_dist, edge_dist * -1))
  
      # plot(int_vims)

  } else {
    
    int_df <- syn_sg_fin %>% 
      dplyr::filter(sampyr == i) %>% 
      mutate(edge_dist = NA) %>% 
      st_join(.,vims_sg_proc %>% dplyr::filter(year == i)) 
  }

  # 
  # plot(int_df %>% dplyr::select(edge_dist))
  # plot(int_vims, add = T)
  
  if (i %in% c(2012, 2014, 2016)){
    int_df %<>% mutate(vims_dens = NA)
  } else if (i == 2013){
    int_df %<>% mutate(vims_dens = ifelse(str_detect(SiteName, "HI"),
                                         NA, 
                                         vims_dens))
  }

  assign('sg_bound', bind_rows(int_df, sg_bound))
}

# process output
sg_bound %<>% 
  mutate(edge_dist = ifelse(abs(edge_dist) > 2500, NA, edge_dist)) %>% 
  dplyr::select(-year, -density) %>% 
  dplyr::rename(density = DENSITY,
                shoots = SHOOTS,
                year = sampyr,
                longitude = Easting__X,
                latitude = Northing__) %>% 
  mutate(year = as.factor(year),
         sg_pres = factor(ifelse(vims_dens > 0 | shoots > 0, 1, 0)),
         vims_pres = factor(ifelse(vims_dens > 0, 1, 0)),
         uva_pres = factor(ifelse(shoots> 0, 1, 0)))

save(sg_bound, file = here::here("data/bound_sg.rdata"))

ggplot(sg_bound) +
  geom_point(aes(x = shoots, y = elevation)) +
  geom_smooth(aes(x = shoots, y = elevation),
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
