library(tidyverse)
library(sf)
library(raster)
library(sdmTMB)
library(patchwork)
library(INLA)
library(ggtext)

# Get seagrass data from presentation-------------------------
load(here::here("data/bound_sg.rdata"))

# Center and rescale function
center.and.rescale <- function(x) {
  x <- x - mean(x, na.rm = T)
  x <- x / (1 * sd(x))
  return(x)
}

sg_df <- sg_bound %>% 
  mutate(meadow = ifelse(str_detect(SiteName, "HI"), "HI", "SB")) %>% 
  dplyr::select(shoots = shoots_raw, 
                elevation = depth, meadow,
                longitude, latitude,
                year, SiteName) %>% 
  st_set_geometry(NULL) %>% 
  filter(year == c(2017), meadow == "HI") %>% 
  mutate(shoots_norm = center.and.rescale(shoots))

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
  ggplot(dem_hi) +
  geom_tile(aes(x = longitude, y = latitude, fill = elevation)) +
  scale_fill_viridis_c()
elevs
# Convert the HI polygon into an INLA mesh------------------------
hi_spat <- hi %>% as("Spatial") %>% inla.sp2segment()

# extract coordinate matrix
coords <- 
  sg_df %>% 
  dplyr::select(longitude, latitude) %>% 
  as.matrix.data.frame()

# Generate the mesh in INLA
hi_mesh <- inla.mesh.2d(coords,
                        boundary=hi_spat, max.edge = 150)

# Last prep step: use sdmTMB to convert INLA mesh into friendly format
hi_spde <- make_mesh(sg_df, c("longitude","latitude"), mesh = hi_mesh)
plot(hi_spde)

# Fit the model----------------------- PROBLEM - model won't converge

# Solution: Bring back the zeros, use a tweedie distribution
# spat_mod <- sdmTMB(shoots ~ elevation,
#                    data = sg_df,
#                    spde = hi_spde,
#                    spatial_only = TRUE)

spat_mod2 <- sdmTMB(shoots ~ elevation,
                    data = sg_df,
                    spde = hi_spde,
                    spatial_only = T)

spat_mod3 <- sdmTMB(shoots ~ s(elevation),
                    data = sg_df,
                    spde = hi_spde,
                    spatial_only = T)

spat_mod4 <- sdmTMB(shoots ~ s(elevation),
                    data = sg_df,
                    spde = hi_spde,
                    spatial_only = T,
                    family = tweedie())

AIC(spat_mod2, spat_mod3, spat_mod4)

# Each basis function gets an entry in the summary
summary(spat_mod4)

# Evaluate the model-----------------
sg_df$res <- residuals(spat_mod4)
qqnorm(sg_df$res);abline(a = 0, b = 1)
hist(sg_df$res)

ggplot(sg_df) +
  geom_point(aes(x = longitude, y= latitude, 
                 size = res, color = res)) +
  labs(size = "Standardized\n residuals",
       color = "Standardized\n residuals") +
  scale_color_continuous(limits=c(-3, 3), breaks=seq(-3, 3, by=1)) +
  scale_size_continuous(limits=c(-3, 3), breaks=seq(-3, 3, by=1)) +
  guides(color= guide_legend(), size=guide_legend()) +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) 

# Predict from the model-------------------
pred_df <- predict(spat_mod4, newdata = dem_hi)

(twed_spat_pred <- 
  ggplot(pred_df) +
  geom_tile(aes(x = longitude, y = latitude, fill = exp(est))) +
  geom_point(data = sg_df, aes(x = longitude, 
                               y = latitude, 
                               size = shoots)) +
  labs(fill = "Predicted seagrass\n density",
       size = "Observed density") +
  scale_fill_viridis_c() +
  theme_bw())


sg_df_pred <- predict(spat_mod4)
(twed_pred <- 
  ggplot(sg_df_pred) +
  geom_point(aes(x = elevation, y= shoots)) +
  geom_line(aes(x = elevation, y = exp(est))) +
  theme_bw() +
  labs(y = "Shoot m<sup>-2</sup>",
       x = "Elevation (m, NAVD88)") +
  theme(axis.title.y = element_markdown(),
        axis.title.x = element_markdown()))

twed_pred + 
  twed_spat_pred 

elevs
hist(dem_hi$elevation)