library(tidyverse)
library(magrittr)
library(sf)
library(raster)
library(sdmTMB)
library(patchwork)
library(INLA)
library(ggtext)
library(fasterize)
library(mapdata)
library(sp)
library(maptools)
library(gstat)
library(rgeos)

# crs
coor_rs <- "+proj=utm +zone=18 +datum=NAD83 +units=km +no_defs"

# Get seagrass data 
load(here::here("data/bound_sg.rdata"))
load(here::here("data/vims_sg.rdata"))

# get VCR polygons
vcr <- st_read(here::here("data/VCR_poly4.kml")) %>% 
  sf::st_transform(crs = st_crs(coor_rs)) %>% 
  st_zm() %>% 
  st_as_sf()

hi <- st_read(here::here("data/HogIslandBay.kml")) %>% 
  sf::st_transform(crs = st_crs(coor_rs)) %>% 
  st_zm() %>% 
  st_as_sf()

# intersect all VIMS imagery with HI polygon
vims_sg_proc %<>% 
  sf::st_transform(crs = st_crs(coor_rs)) %>% 
  st_intersection(.,hi)

# read in DEM
dem_raw <- raster(here::here("data/V2_topobathy_tif.tif"))
# dem_raw <- projectRaster(dem_raw, crs = coor_rs)
dem <- mask( crop( dem_raw,y = extent(hi %>% 
                                        st_transform("+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs")) ), hi %>% 
               st_transform("+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs"))
dem[dem > 0 ] <- NA
dem <- aggregate(dem, fact = 8.5, fun = mean)
dem <- projectRaster(dem, crs = coor_rs)

# Center and rescale function
center.and.rescale <- function(x) {
  x <- x - mean(x, na.rm = T)
  x <- x / (1 * sd(x))
  return(x)
}

# prediction function
plot_map <- function(dat, column) {
  ggplot(dat, aes_string("x", "y", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}

# process seagrass data (and filter to HI only)
sg_df <- sg_bound %>% 
  filter(meadow == "HI") %>% 
  st_transform(st_crs(coor_rs)) %>% 
  {. ->> sg_sf} %>% 
  dream::sfc_as_cols() %>% 
  st_set_geometry(NULL) %>% 
  dplyr::rename(shoots = SHOOTS) %>% 
  mutate(shoots_norm = center.and.rescale(shoots),
         year = as.numeric(as.character(year)),
         bathy = raster::extract(dem, y = sg_sf %>% as_Spatial()))


rast_res <- 0.025
# template raster
r <- raster(extent(hi),
            res = rast_res,
            crs =coor_rs)

# handle seagrass aerial imagery
vcr_sg_df <- NULL
for (i in unique(vims_sg_proc$year)){
  if (i %in% c(2014, 2016)){
    missing <- 
      rasty %>% 
      as("SpatialPixelsDataFrame") %>% 
      as_tibble() %>% 
      dplyr::rename(vims_dens = layer) %>% 
      mutate(year = as.numeric(i))
    assign("vcr_sg_df",
           rbind(missing, 
                 vcr_sg_df))
    next
  } else if (i == 2013){
    int_df %<>% 
      mutate(vims_dens = NA,
             year = as.numeric(i))
  } else {
    rasty <- 
      vims_sg_proc %>% 
      filter(year == i, vims_dens > 0) %>% 
      st_cast("MULTIPOLYGON") %>%
      fasterize::fasterize(raster = r, background = 0, field = "vims_dens")
    
    int_df <- rasty %>% 
      as("SpatialPixelsDataFrame") %>% 
      as_tibble() %>% 
      mutate(year = as.numeric(i)) %>% 
      dplyr::rename(vims_dens = layer) 
  }
  
  assign("vcr_sg_df",
         rbind(int_df,
               vcr_sg_df))
  
}

vcr_sg_df$bathy <- rep(raster::extract(dem, rasty %>% 
                  as("SpatialPixelsDataFrame")),
      length(unique(vims_sg_proc$year)))

# ggplot(vcr_sg_df) +
#   geom_tile(aes(y = y, x = x, fill = vims_dens)) +
#   facet_wrap(~year)

# read polygons for mesh creation
file_vec <- list.files(here::here("data/vcr_barrier_model_polygons/"))
bar_poly <- NULL
for (i in file_vec){
  assign("bar_poly", 
         bind_rows(st_read(here::here("data/vcr_barrier_model_polygons",i)) %>% 
                     sf::st_transform(crs = st_crs(coor_rs)) %>% 
                     st_zm() %>% 
                     st_as_sf() %>% 
                     mutate(water = ifelse(i == "vcr_mesh_outer1.kml", 1, 0)),
                   bar_poly))
}

# create barrier mesh
land_sp <- 
  bar_poly %>% 
  filter(water == 0) %>% 
  summarise(do_union = F)  %>% 
  st_intersection(.,hi) %>% 
  as("Spatial")
land_sf <- 
  bar_poly %>% 
  filter(water == 0) %>% 
  summarise(do_union = T) %>% 
  st_intersection(.,hi)
vcr_sp <- 
  vims_sg_proc %>% 
  filter(year == 2018, vims_dens > 0) %>% 
  dplyr::select(vims_dens) %>% 
  summarise() %>% 
  st_buffer(.,0.55) %>% 
  st_difference(land_sf) %>% 
  {. ->> bar_sf} %>% 
  as("Spatial")
(bar_poly2 <- rgeos::gDifference(vcr_sp, land_sp))
bar_rast <- fasterize(bar_sf, raster = r)

ggplot() +
  geom_point(data = sg_df %>% filter(year == 2007), aes(x = x, y = y)) +
  geom_sf(data = bar_sf)

max.edge = 0.25
bound.outer = 0.5
mesh <- inla.mesh.2d(boundary = bar_poly2,
                     max.edge = c(1,5) * max.edge,
                     cutoff = 0.05,
                     offset = c(max.edge, bound.outer))

water.tri = inla.over_sp_mesh(bar_poly2, y = mesh, 
                              type = "centroid", ignore.CRS = TRUE)
num.tri = length(mesh$graph$tv[, 1])
barrier.tri = base::setdiff(1:num.tri, water.tri)
poly.barrier = inla.barrier.polygon(mesh, barrier.triangles = barrier.tri)
plot(poly.barrier)


yrs <- 2007:2012

mod_df <-  sg_df %>% filter(!is.na(x),
                            year %in% yrs) %>% 
  mutate(vims_pres = as.numeric(as.character(vims_pres)),
         year = as.numeric(as.character(year)))

vcr_spde <- sdmTMB::make_mesh(data = mod_df, xy_cols = c("x","y"),mesh = mesh)
# vcr_spde <- add_barrier_mesh(vcr_spde, land_sf)

# fit model
mod <- 
  sdmTMB(shoots ~ 0 + vims_dens + s(bathy, k = 4) + factor(year), 
         data = mod_df, 
         spde = vcr_spde,
         # extra_time = c(2013L, 2014L, 2016L),
         # fields = "AR1",
         time = "year",
         family = tweedie(link = "log"))
summary(mod)
hist(residuals(mod))


t <- raster::extract(bar_rast, vcr_sg_df %>% 
                       filter(year == 2007) %>% 
                       .[c("x","y")])
pred_df <- vcr_sg_df %>% 
  filter(year %in% yrs) %>% 
  mutate(pres = rep(t, length(yrs)))

predictions <- predict(mod, newdata = pred_df)
plot_map(predictions %>% 
           filter(!is.na(pres),
                  vims_dens != 0), "exp(est)") +
  # scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)") 

nd <- tibble(bathy = seq(min(mod_df$bathy), max(mod_df$bathy), 
                         length.out = 100),
         vims_dens = mean(mod_df$vims_dens),
         year = 2012)

sg_df_pred <- predict(mod, newdata = nd, type = "response", se_fit = TRUE, re_form = NA)
(twed_pred <- 
    ggplot(sg_df_pred) +
    geom_point(data = mod_df, aes(x = bathy, y= shoots)) +
    geom_line(aes(x = bathy, y = exp(est))) +
    theme_bw() +
    labs(y = "Shoot canopy height (cm)",
         x = "Elevation (m, NAVD88)") +
    theme(axis.title.y = element_markdown(),
          axis.title.x = element_markdown()))

nd <- tibble(vims_dens = seq(min(mod_df$vims_dens), max(mod_df$vims_dens), 
                         length.out = 100),
             bathy = mean(mod_df$bathy),
             year = 2012)

sg_df_pred <- predict(mod, newdata = nd, type = "response", se_fit = TRUE, re_form = NA)
(twed_pred <- 
    ggplot(sg_df_pred) +
    geom_point(data = mod_df, aes(x = vims_dens, y= shoots)) +
    geom_line(aes(x = vims_dens, y = exp(est))) +
    theme_bw() +
    labs(y = "Shoot canopy height (cm)",
         x = "Elevation (m, NAVD88)") +
    theme(axis.title.y = element_markdown(),
          axis.title.x = element_markdown()))

pdf <- predict(mod)
ggplot(pdf) +
  geom_point(aes(y = exp(est), x = shoots)) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~year)
