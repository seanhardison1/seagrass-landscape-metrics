library(tidyverse)
library(sf)
library(raster)
library(sdmTMB)
library(patchwork)
library(INLA)
library(ggtext)
library(fasterize)

# Get seagrass data from presentation-------------------------
load(here::here("data/bound_sg.rdata"))
load(here::here("data/HogIslandBaySeagrass.rdata"))
load(here::here("data/processed_dem.rdata"))

# Center and rescale function
center.and.rescale <- function(x) {
  x <- x - mean(x, na.rm = T)
  x <- x / (1 * sd(x))
  return(x)
}

sg_df <- sg_bound %>% 
  mutate(meadow = ifelse(str_detect(SiteName, "HI"), "HI", "SB")) %>% 
  st_set_geometry(NULL) %>% 
  filter(meadow == "HI") %>% 
  mutate(shoots_norm = center.and.rescale(shoots),
         year = as.numeric(as.character(year)),
         bpm = shoots * biomass) %>% 
  filter(!is.na(biomass))

r <- raster(extent(hi_sg),
            res = 1,
            crs ="+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")
hi_sg_df <- NULL
area_df <- NULL
for (i in unique(hi_sg$YEAR)){
  rasty <- 
    hi_sg %>% 
      filter(YEAR == i) %>% 
      fasterize::fasterize(raster = r, background = 0)
  
  int_df <- rasty %>% 
    as("SpatialPixelsDataFrame") %>% 
    as_tibble() %>% 
    mutate(year = i) %>% 
    dplyr::rename(longitude = x, latitude = y) %>% 
    mutate(elevation = raster::extract(dem_resample, 
                                     .[c("longitude","latitude")]))
  assign("hi_sg_df",
         rbind(int_df,
               hi_sg_df))

}

area_df <-  hi_sg %>% 
              dplyr::group_by(YEAR) %>% 
              summarise(do_union = F) %>% 
              mutate(sg_area =  as.numeric(st_area(.))) %>% 
              st_set_geometry(NULL)


# Read HI data------------------------
hi <- st_read(here::here("data/HI_INLA_poly.kml")) %>% 
  sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")) %>% 
  st_zm()


# Read and crop raster------------------------
load(here::here("data/processed_dem.rdata"))
dem <- mask( crop( dem_resample,y = extent(dem_resample) ), hi)

sg_df$elevation <- 
  raster::extract(dem, y = sg_df[,c("longitude","latitude")])

# Get elevation data to predict over from raster------------------------
dem_hi <- dem %>% 
  as("SpatialPixelsDataFrame") %>%
  as_tibble() %>% 
  dplyr::rename(elevation = 1,
                longitude = 2,
                latitude = 3)

elevs <- 
  ggplot() +
  geom_sf(data = hi) +
  geom_tile(data = hi_sg_df, aes(x = longitude, y = latitude)) +
  # geom_tile(data = dem_hi,
  #           aes(x = longitude, y = latitude, fill = elevation)) +
  scale_fill_viridis_c()
# elevs

# Convert the HI polygon into an INLA mesh------------------------
hi_spat <- hi %>% as("Spatial") %>% inla.sp2segment()

# extract coordinate matrix
coords <- 
  sg_df %>% 
  dplyr::select(longitude, latitude) %>% 
  as.matrix.data.frame()

# Generate the mesh in INLA
hi_mesh <- inla.mesh.2d(coords,
                        boundary=hi_spat, max.edge = 100)

# Last prep step: use sdmTMB to convert INLA mesh into friendly format
hi_spde <- make_mesh(sg_df, c("longitude","latitude"), 
                      n_knots = 35, type = "kmeans")
plot(hi_spde)

# Fit the model----------------------- 
spat_mod4 <- sdmTMB(bpm ~ 0 + as.factor(year),
                    data = sg_df,
                    spde = hi_spde,
                    time = "year",
                    family = tweedie())

spat_mod5 <- sdmTMB(bpm ~ 0 + as.factor(year),
                    data = sg_df,
                    spde = hi_spde,
                    time = "year",
                    ar1_fields = TRUE,
                    nlminb_loops = 2,
                    family = tweedie())

spat_mod6 <- sdmTMB(bpm ~ 0 + as.factor(year),
                    data = sg_df,
                    spde = hi_spde,
                    time = "year",
                    nlminb_loops = 2,
                    family = tweedie())

AIC(spat_mod4, spat_mod5, spat_mod6)

# Each basis function gets an entry in the summary
summary(spat_mod6)

# Evaluate the model-----------------
# sg_df$res <- residuals(spat_mod6)
qqnorm(residuals(spat_mod6));abline(a = 0, b = 1)
# hist(sg_df$res)


# Predict from the model-------------------
ndata <- hi_sg_df %>% 
  filter(year %in% unique(sg_df$year))

pred_df <- predict(spat_mod6, newdata = ndata) %>% 
  mutate_at(vars(est:epsilon_st), function(x){.$layer * x}) %>% 
  filter(est > 0)

(twed_spat_pred <- 
    ggplot(pred_df) +
    geom_tile(aes(x = longitude, y = latitude, fill = exp(est) * 100)) +
    labs(fill = "Predicted seagrass\n biomass (g m<sup>-2</sup>)") +
    scale_fill_viridis_c() +
    facet_wrap(~year) +
    theme_bw() +
    theme(legend.title = element_textbox(),
          legend.position = "bottom",
          legend.text = element_text(angle = 45)))

pred_df %>% 
  mutate(est = exp(est) * 100) %>% 
  group_by(year) %>% 
  summarise(est = sum(est)) %>% 
  inner_join(.,area_df) %>% 
  mutate(biomass = (area * est))

sg_df_pred <- predict(spat_mod6)
(twed_pred <- 
    ggplot(sg_df_pred) +
    geom_point(aes(x = elevation, y= bpm)) +
    geom_line(aes(x = elevation, y = exp(est))) +
    theme_bw() +
    labs(y = "Shoot m<sup>-2</sup>",
         x = "Elevation (m, NAVD88)") +
    theme(axis.title.y = element_markdown(),
          axis.title.x = element_markdown()) +
    facet_wrap(~year))
