library(tidyverse)
library(magrittr)
library(sf)
library(raster)
library(sdmTMB)
library(patchwork)
library(INLA)
library(ggtext)
library(fasterize)
library(INLA)
library(mapdata)
library(sp)
library(maptools)
library(gstat)
library(rgeos)

# Get seagrass data from presentation-------------------------
load(here::here("data/bound_sg.rdata"))
load(here::here("data/vims_sg.rdata"))
vims_sg_proc %<>% sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=NAD83 +units=km +no_defs"))


vcr <- st_read(here::here("data/VCR_poly4.kml")) %>% 
  sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=NAD83 +units=km +no_defs")) %>% 
  st_zm() %>% 
  st_as_sf()

# Center and rescale function
center.and.rescale <- function(x) {
  x <- x - mean(x, na.rm = T)
  x <- x / (1 * sd(x))
  return(x)
}

sg_df <- sg_bound %>% 
  st_transform(st_crs("+proj=utm +zone=18 +datum=WGS84 +units=km +no_defs")) %>% 
  dream::sfc_as_cols() %>% 
  st_set_geometry(NULL) %>% 
  dplyr::rename(shoots = SHOOTS) %>% 
  mutate(shoots_norm = center.and.rescale(shoots),
         year = as.numeric(as.character(year)))

r <- raster(extent(vcr),
            res = 0.5,
            crs ="+proj=utm +zone=18 +datum=WGS84 +units=km +no_defs")

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
  }
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
  
  if (i == 2013){
    int_df %<>% 
      mutate(vims_dens = ifelse(y > 4127.5, NA, vims_dens))
  }
  
  assign("vcr_sg_df",
         rbind(int_df,
               vcr_sg_df))

}

ggplot(vcr_sg_df) +
  geom_tile(aes(y = y, x = x, fill = vims_dens)) +
  facet_wrap(~year)

# read polygons for mesh creation
file_vec <- list.files(here::here("data/vcr_barrier_model_polygons/"))
bar_poly <- NULL
for (i in file_vec){
  assign("bar_poly", 
         bind_rows(st_read(here::here("data/vcr_barrier_model_polygons",i)) %>% 
                     sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=NAD83 +units=km +no_defs")) %>% 
                     st_zm() %>% 
                     st_as_sf() %>% 
                     mutate(water = ifelse(i == "vcr_mesh_outer1.kml", 1, 0)),
                   bar_poly))
}

# create barrier mesh
land_sp <- 
  bar_poly %>% 
  filter(water == 0) %>% 
  summarise(do_union = F) %>% 
  as("Spatial")
land_sf <- 
  bar_poly %>% 
  filter(water == 0) %>% 
  summarise(do_union = F)
vcr_sp <- 
  bar_poly %>% 
  filter(water == 1) %>% 
  summarise(do_union = F) %>% 
  as("Spatial")
bar_poly2 <- rgeos::gDifference(vcr_sp, land_sp)

max.edge = 1
bound.outer = 0.5
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


yrs <- 2007:2012

mod_df <-  sg_df %>% filter(!is.na(x),
                            year %in% yrs) %>% 
  mutate(vims_pres = as.numeric(as.character(vims_pres)))

vcr_spde <- sdmTMB::make_mesh(data = mod_df, xy_cols = c("x","y"),mesh = mesh)
vcr_spde <- add_barrier_mesh(vcr_spde, land_sf)

# sdm
mod <- 
  sdmTMB(shoots ~ vims_dens, 
         data = mod_df, 
         spde = vcr_spde,
         # extra_time = c(2013L, 2014L, 2016L),
         # fields = "AR1",
         time = "year",
         family = tweedie(link = "log"))
summary(mod)
hist(residuals(mod))

# 
# p2013 <- vcr_sg_df %>% filter((year == 2014)) %>% mutate(year = 2013)
# p_2014_2015 <- vcr_sg_df %>% filter(year %in% c(2014, 2016))
# interp_pred_df <- rbind(p2013, p_2014_2015)


plot_map(predictions, "exp(est)") +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")

predictions <- predict(mod, newdata = vcr_sg_df)

complete <- predictions %>% 
  bind_rows(., 
            interp_predictions)


# hurdle model
m_bin <- sdmTMB(vims_pres ~ 1,
                fields = "AR1",
                data = mod_df, time = "year", spde = vcr_spde,
                family = binomial(link = "logit"), include_spatial = FALSE)
summary(m_bin)
hist(residuals(m_bin))

p2013 <- vcr_sg_df %>% filter(year %in% 2007:2009)

interp_predictions <- predict(m_bin, newdata = p2013)
plot_map(interp_predictions %>% filter(vims_dens > 0 ), "plogis(est)") +
  scale_fill_viridis_c()
