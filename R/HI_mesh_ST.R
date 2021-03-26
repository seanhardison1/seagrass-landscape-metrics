library(tidyverse)
library(sf)
library(raster)
library(sdmTMB)
library(patchwork)
library(INLA)
library(ggtext)

# Get seagrass data from presentation-------------------------
load(here::here("data/bound_sg.rdata"))


sg_df <- sg_bound %>% 
  mutate(meadow = ifelse(str_detect(SiteName, "HI"),
                         "HI", "SB")) %>% 
  # dplyr::select(elevation = depth, meadow,
  #               longitude, latitude,
  #               year, SiteName) %>% 
  st_set_geometry(NULL) %>% 
  filter(meadow == "HI", year == 2017) 


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


sg_df2 <- sg_df %>% filter(!is.na(canopy))

# extract coordinate matrix
coords <- 
  sg_df2 %>% 
  dplyr::select(longitude, latitude) %>% 
  as.matrix.data.frame()

# Generate the mesh in INLA
hi_mesh <- inla.mesh.2d(coords, 
                        boundary=hi_spat, max.edge = 150)
plot(hi_mesh)
# Last prep step: use sdmTMB to convert INLA mesh into friendly format
hi_spde <- make_mesh(sg_df2, c("longitude","latitude"),
                     mesh = hi_mesh)
plot(hi_spde)

# Fit the model----------------------- PROBLEM - model won't converge
ggplot(sg_df2) +
  geom_point(aes(x = elevation, y= canopy))

# Solution: Bring back the zeros, use a tweedie distribution
spat_mod <- sdmTMB(canopy ~ elevation,
                   data = sg_df2,
                   spde = hi_spde,
                   spatial_only = TRUE)

spat_mod2 <- sdmTMB(canopy ~ s(elevation),
                    data = sg_df2,
                    spde = hi_spde,
                    spatial_only = T)

spat_mod3 <- sdmTMB(shoots ~ s(canopy),
                    data = sg_df,
                    spde = hi_spde,
                    spatial_only = T,
                    family = tweedie())

AIC(spat_mod, spat_mod2)

# Each basis function gets an entry in the summary
summary(spat_mod)

# Evaluate the model-----------------

sg_df2$res <- residuals(spat_mod)
qqnorm(sg_df2$res);abline(a = 0, b = 1)
hist(sg_df2$res)

ggplot(sg_df2) +
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
pred_df <- predict(spat_mod2, 
                   newdata = dem_hi)

(twed_spat_pred <- 
  ggplot(pred_df) +
  geom_tile(aes(x = longitude, y = latitude, fill = est)) +
  geom_point(data = sg_df2, aes(x = longitude, 
                               y = latitude, 
                               size = canopy)) +
  labs(fill = "Predicted seagrass\n canopy height (cm)",
       size = "Observed canopy\n height (cm)") +
  scale_fill_viridis_c() +
  theme_bw())


sg_df_pred <- predict(spat_mod)
(twed_pred <- 
  ggplot(sg_df_pred) +
  geom_point(aes(x = elevation, y= canopy)) +
  geom_line(aes(x = elevation, y = est)) +
  theme_bw() +
  labs(y = "Shoot canopy height (cm)",
       x = "Elevation (m, NAVD88)") +
  theme(axis.title.y = element_markdown(),
        axis.title.x = element_markdown()))

twed_pred + 
  twed_spat_pred 

 elevs
hist(dem_hi$elevation)
