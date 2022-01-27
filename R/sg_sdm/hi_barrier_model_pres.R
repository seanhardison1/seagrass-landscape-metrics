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

# for visualizing
vcr_yr <- 
  vims_sg_proc %>% 
  filter(vims_dens > 0) %>% 
  dplyr::select(vims_dens, year) %>% 
  group_by(year) %>% 
  summarise()

# read in DEM
dem_raw <- raster(here::here("data/V2_topobathy_tif.tif"))
# dem_raw <- projectRaster(dem_raw, crs = coor_rs)
dem <- mask( crop( dem_raw,y = extent(hi %>% 
                                        st_transform("+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs")) ), hi %>% 
               st_transform("+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs"))
dem[dem > 0 ] <- NA
dem <- aggregate(dem, fact = 8.5, fun = mean)
dem <- projectRaster(dem, crs = coor_rs)

# prediction function
plot_map <- function(dat, column, vcr_sp = vcr_yr) {
  ggplot() +
    geom_raster(data = dat, aes_string("x", "y", fill = column)) +
    geom_sf(data = vcr_sp, fill = "transparent", color = "purple") +
    facet_wrap(~year) +
    coord_sf()  +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "bottom") +
    scale_fill_gradientn(colors = pals::ocean.deep(10)) +
    ecodata::theme_map()
    
}

rast_res <- 0.1
# template raster
r <- raster(extent(hi),
            res = rast_res,
            crs =coor_rs)


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
  st_buffer(.,0.1) %>% 
  st_difference(land_sf) %>% 
  {. ->> bar_sf} %>% 
  as("Spatial")
(bar_poly2 <- rgeos::gDifference(vcr_sp, land_sp))
bar_rast <- fasterize(bar_sf, raster = r)

max.edge = 0.1
bound.outer = 0.75
mesh <- inla.mesh.2d(boundary = bar_poly2,
                     max.edge = c(1,5) * max.edge,
                     cutoff = 0.05,
                     offset = c(max.edge, bound.outer))
plot(mesh)
water.tri = inla.over_sp_mesh(bar_poly2, y = mesh, 
                              type = "centroid", ignore.CRS = TRUE)
num.tri = length(mesh$graph$tv[, 1])
barrier.tri = base::setdiff(1:num.tri, water.tri)
poly.barrier = inla.barrier.polygon(mesh, barrier.triangles = barrier.tri)
plot(poly.barrier)


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
      fasterize::fasterize(raster = r, background = 0, field = "vims_dens") %>% 
      mask(.,bar_rast) 
      
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


ggplot(vcr_sg_df) +
  geom_raster(aes(y = y, x = x, fill = vims_dens)) +
  facet_wrap(~year)



yrs <- c(2017)
# yrs <- 2007:2012

mod_df <-  vcr_sg_df %>% filter(year %in% yrs) %>% 
  mutate(vims_pres = ifelse(vims_dens > 0 ,1, 0))

vcr_spde <- sdmTMB::make_mesh(data = mod_df, xy_cols = c("x","y"),mesh = mesh)

plot(vcr_spde)

mod <-
  sdmTMB(vims_pres ~ 1,
         data = mod_df,
         spde = vcr_spde,
         # extra_time = c(2013L, 2014L, 2016L),
         # fields = "AR1",
         # previous_fit = mod,
         # time = "year",
         # spatial_trend = TRUE,
         family = binomial())


# mod <-
#   sdmTMB(vims_pres ~ 1,
#          data = mod_df,
#          spde = vcr_spde,
#          extra_time = c(2013L, 2014L, 2016L),
#          fields = "AR1",
#          # previous_fit = mod,
#          time = "year",
#          # spatial_trend = TRUE,
#          family = binomial())

vcr_spde_bar <- add_barrier_mesh(vcr_spde, bar_sf)
bar_mod <- 
  sdmTMB(vims_pres ~ 1, 
         data = mod_df, 
         spde = vcr_spde_bar,
         extra_time = c(2013L, 2014L, 2016L),
         fields = "AR1",
         # previous_fit = mod,
         time = "year",
         # spatial_trend = TRUE,
         family = binomial())

AIC(mod, bar_mod)
summary(mod)
summary(bar_mod)
hist(residuals(mod))
t <- raster::extract(bar_rast, vcr_sg_df %>% 
                       filter(year == 2007) %>% 
                       .[c("x","y")])
pred_df <- vcr_sg_df %>% 
  filter(year %in% yrs) %>% 
  mutate(pres = rep(t, length(yrs)))

predictions <- predict(mod, newdata = pred_df)
plot_map(dat = predictions, column =  "plogis(est)") +
  # scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)") 

bar_predictions <- predict(bar_mod, newdata = pred_df)
plot_map(dat = bar_predictions, column =  "plogis(est)") +
  # scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Barrier Prediction (fixed effects + all random effects)") 

ggplot() +
  geom_raster(data = vcr_sg_df%>% 
                filter(year %in% yrs), aes(x = x, y = y, fill = vims_dens), color = "red") +
  facet_wrap(~year)

fc_pred_df <- pred_df %>% 
  bind_rows(pred_df %>% filter(year == 2012) %>% mutate(year = 2013),
            pred_df %>% filter(year == 2012) %>% mutate(year = 2014),
            pred_df %>% filter(year == 2012) %>% mutate(year = 2016)) %>% 
  bind_rows(vcr_sg_df%>% 
    filter(year != 2016))
unique(fc_pred_df$year)

bar_mod_fc <- predict(bar_mod, newdata = fc_pred_df)

spat_int_fig <- 
  plot_map(dat = bar_mod_fc %>% 
           mutate(year = ifelse(year %in% c(2013, 2014, 2016),
                                paste(year, "(int.)"),
                                year)), column =  "plogis(est)") +
  # scale_fill_viridis_c(trans = "sqrt") +
  labs(fill = "Eelgrass\n occurrence prob.",
       title = "Spatiotemporal interpolation of eelgrass presence",
       subtitle = "Hog Island Bay (2007-2018)")

ggsave(spat_int_fig, filename = here::here("figures/spat_int_fig.png"),
       dpi = 200, width = 6, height = 8)


mod_fc <- predict(mod, newdata = fc_pred_df)
plot_map(dat = mod_fc, column =  "plogis(est)") +
  # scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Non-barrier model forecast (fixed effects + all random effects)") 

bar_mod_fc_bin <- bar_mod_fc
mod_fc_bin <- mod_fc
save(bar_mod_fc_bin, mod_fc_bin, file = here::here("data/stage1_predictions.rdata"))
