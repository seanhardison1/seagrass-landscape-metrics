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

#load predictions from binomial model
load(here::here("data/stage1_predictions.rdata"))

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
# 
ggplot() +
  geom_sf(data = vcr_yr %>%
            filter(year %in% yrs),
          fill = "green", alpha = 0.5) +
  geom_sf(data =barrier_sf, alpha = 0.5, 
          fill = "purple") +
  facet_wrap(~year) +
  geom_point(data = mod_df_dens %>%
               filter(shoots > 10,
                      year %in% yrs),
             aes(x = x, y = y,
                 color = shoots),
             size = 0.5,
             color = "black") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank()) +
  scale_color_distiller(palette = "RdBu",
                        type = "div") 




rast_res <- 0.01
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
  st_buffer(.,0.05) %>% 
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
barrier_sf <- as(poly.barrier, "sf")
st_crs(barrier_sf) <- coor_rs

mod_df_dens <- 
  sg_bound %>% 
  st_transform(crs = "+proj=utm +zone=18 +datum=NAD83 +units=km +no_defs") %>% 
  dream::sfc_as_cols() %>% 
  st_set_geometry(NULL) %>% 
  filter(meadow == "HI") %>% 
  dplyr::rename(site = SiteName,
                shoots = SHOOTS) %>% 
  mutate(year = as.numeric(as.character(year)))

mod_df_dens %>% 
  group_by(year) %>% 
  dplyr::summarise(n = n())

ggplot(mod_df_dens) +
  geom_point(aes(x = x, y = y, color = shoots)) +
  facet_wrap(~year)

yrs <- c(2017, 2018)
mod_df_dens2 <- mod_df_dens %>% filter(year %in% yrs)
vcr_spde <- sdmTMB::make_mesh(data =mod_df_dens2, xy_cols = c("x","y"),mesh = mesh)
vcr_spde_bar <- add_barrier_mesh(vcr_spde, bar_sf)
plot(vcr_spde_bar)

# AR1

bar_mod_dens_ar1 <- 
  sdmTMB(shoots ~ 
           as.factor(vims_dens), 
         data = mod_df_dens2, 
         spde = vcr_spde_bar,
         spatial = "on",
         time = "year",
         spatiotemporal = "AR1",
         family = tweedie())
# summary(bar_mod_dens)

mod_dens_ar1 <- 
  sdmTMB(shoots ~ 
           as.factor(vims_dens), 
         data = mod_df_dens2, 
         spde = vcr_spde,
         spatial = "on",
         time = "year",
         spatiotemporal = "AR1",
         family = tweedie())

# IID
bar_mod_dens <- 
  sdmTMB(shoots ~ 
           as.factor(vims_dens), 
         data = mod_df_dens2, 
         spde = vcr_spde_bar,
         spatial = "on",
         time = "year",
         spatiotemporal = "IID",
         family = tweedie())
# summary(bar_mod_dens)

mod_dens <- 
  sdmTMB(shoots ~ 
           as.factor(vims_dens), 
         data = mod_df_dens2, 
         spde = vcr_spde,
         spatial = "on",
         time = "year",
         spatiotemporal = "IID",
         family = tweedie())

# model compare
AIC(bar_mod_dens, mod_dens)

bar_rast_df <- 
  bar_rast %>% 
  as("SpatialPixelsDataFrame") %>% 
  as.data.frame() %>% 
  tidyr::expand_grid(vims_dens = unique(mod_df_dens2$vims_dens),
                     year = yrs) %>% 
  dplyr::select(-layer)

bar_mod_fc <- sdmTMB:::predict.sdmTMB(bar_mod_dens, newdata = bar_rast_df)
mod_fc <- sdmTMB:::predict.sdmTMB(mod_dens, newdata = bar_rast_df)

# prediction function
plot_map <- function(dat, column, vcr_sp = vcr_yr) {
  ggplot() +
    geom_raster(data = dat, aes_string("x", "y", fill = column)) +
    geom_sf(data = vcr_sp %>% 
              filter(year %in% unique(dat$year)), fill = "transparent", color = "purple") +
    coord_sf()  +
    facet_wrap(~year) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "bottom") +
    scale_fill_gradientn(colors = pals::ocean.deep(10),
                         breaks = c(0, 200, 400, 600),
                         limits = c(0,600)) +
    theme_minimal()
  
}

plot_map(dat = bar_mod_fc, column = "exp(est)") +
  labs(fill = "Density\n(shoots 10 m<sup>-2</sup>)") +
  theme(legend.title = element_markdown(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")

plot_map(dat = mod_fc, column = "exp(est)") +
  labs(fill = "Density\n(shoots 10 m<sup>-2</sup>)") +
  theme(legend.title = element_markdown(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")
  
AIC(bar_mod_dens, bar_mod_dens_ar1,
    mod_dens, mod_dens_ar1)

# dens_index <- sdmTMB::get_index(bar_mod_fc, bias_correct = TRUE)

dens_index %>% 
  mutate(year = as.numeric(as.character(year))) %>% 
  ggplot() + 
    geom_line(aes(y = est, x = year))

combined_preds <- 
  bar_mod_fc$data %>% 
  mutate(year = as.numeric(as.character(year))) %>% 
  dplyr::rename(dens_est = est) %>% 
  left_join(.,mod_fc$data %>% 
              mutate(year = as.numeric(as.character(year))) %>% 
              dplyr::select(x, y, year, dens_est_sgf = est)) %>% 
  left_join(.,bar_mod_fc_bin %>% 
              dplyr::select(x, y, year, bin_est = est)) %>% 
  left_join(.,mod_fc_bin %>% 
              dplyr::select(x, y, year, bin_est_sgf = est)) %>% 
  mutate(combined_est_bar = plogis(bin_est) * exp(dens_est),
         combined_est_sgf = plogis(bin_est_sgf) * exp(dens_est_sgf))

plot_map(dat = combined_preds, column = "combined_est_bar") +
  labs(fill = "Shoot density",
       title = "Hurdle-barrier model predictions")

plot_map(dat = combined_preds, column = "combined_est_sgf") +
  labs(fill = "Shoot density",
       title = "Hurdle-SGF model predictions")

combined_preds %>% 
  group_by(year) %>% 
  dplyr::summarise(est_bar = sum(combined_est_bar),
                   est_sgf = sum(combined_est_sgf)) %>% 
  ggplot() + 
    geom_line(aes(x = year, y = est_bar), color = "blue") + 
    geom_line(aes(x = year, y = est_sgf))

mod_df_dens2 %>% 
  mutate(year = as.numeric(as.character(year))) %>% 
  group_by(year) %>% 
  dplyr::summarise(shoots = mean(shoots)) %>% 
  ggplot() + 
    geom_line(aes(x = year, y = shoots))

