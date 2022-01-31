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
library(mgcv)
load(here::here("data/iid_mod_output.rdata"))
# prediction function
plot_map <- function(dat, column, vcr_sp = vcr_yr, sdmtmb = T) {
  if(sdmtmb){
    msd <- max(exp(dat$est), na.rm = T)
    } else{
    msd <- max(dat$fit, na.rm = T)
    }
  
  ggplot() +
    geom_raster(data = dat, aes_string("x", "y", fill = column)) +
    geom_sf(data = vcr_sp %>% 
              filter(year %in% unique(dat$year)), fill = "transparent", color = "purple") +
    coord_sf()  +
    # facet_wrap(~year) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "bottom") +
    scale_fill_gradientn(colors = pals::ocean.deep(10),
                         limits = c(0,msd)) +
    theme_minimal()
  
}

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

yrs <- c(2007, 2008, 2009, 2010, 2011, 2012, 2015, 2017, 2018)
pred_shoot_dens <- NULL
shoot_index <- NULL
plt_list <- list()
mod_list <- list()
for (yr in 1:length(yrs)){
  print(yrs[yr])
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
    filter(year == yrs[yr], 
           vims_dens > 0) %>% 
    dplyr::select(vims_dens) %>%
    {. ->> pred_poly} %>% 
    summarise() %>% 
    st_buffer(.,0.025) %>% 
    st_difference(land_sf) %>% 
    {. ->> bar_sf} %>% 
    as("Spatial")
  bar_poly2 <- rgeos::gDifference(vcr_sp, land_sp)

  if (!yr %in% c(8,9)){
    mesh <- inla.mesh.2d(boundary = bar_poly2,
                         max.edge = c(0.02, 1),
                         cutoff = 0.05,
                         offset = c(0.05, 0.5))
  } else {
    mesh <- inla.mesh.2d(boundary = bar_poly2,
                         max.edge = c(0.05, 1),
                         cutoff = 0.05,
                         offset = c(0.05, 0.5))
  }

  plot(mesh)
  water.tri = inla.over_sp_mesh(bar_poly2, y = mesh, 
                                type = "centroid", ignore.CRS = TRUE)
  num.tri = length(mesh$graph$tv[, 1])
  barrier.tri = base::setdiff(1:num.tri, water.tri)
  poly.barrier = inla.barrier.polygon(mesh, barrier.triangles = barrier.tri)
  
  # plot(poly.barrier)
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
    mutate(year = dream::fct_to_num(year))
  
  mod_df_dens2 <- mod_df_dens[mod_df_dens$year == yrs[yr],]
  vcr_spde <- sdmTMB::make_mesh(data =mod_df_dens2, xy_cols = c("x","y"),mesh = mesh)
  vcr_spde_bar <- add_barrier_mesh(vcr_spde, bar_sf)
  # plot(vcr_spde_bar)
  
  # Barrier
  bar_mod_dens <- 
    sdmTMB(shoots ~ 
             vims_dens, 
           data = mod_df_dens2, 
           mesh = vcr_spde_bar,
           spatial = "on",
           reml = TRUE,
           family = tweedie())
  
  # Stationary
  mod_dens <- 
    sdmTMB(shoots ~ 
             vims_dens, 
           data = mod_df_dens2, 
           mesh = vcr_spde,
           spatial = "on",
           reml = TRUE,
           family = tweedie())
    
  # No boundary considerations
  if (yrs[yr] %in% 2007){
    mod_dens_gam <- 
      gam(shoots ~ vims_dens + s(x, y, bs = "gp", k = 4),
          family = "tw",
          method = "REML",
          data = mod_df_dens2 )
  } else {
    mod_dens_gam <- 
      gam(shoots ~ vims_dens + s(x, y, bs = "gp"),
          family = "tw",
          method = "REML",
          data = mod_df_dens2 )
  }
  
  
  AIC(bar_mod_dens, mod_dens, mod_dens_gam)
  
  bar_rast <- fasterize(pred_poly, raster = r, field = "vims_dens")
  bar_rast_df <- 
    bar_rast %>% 
    as("SpatialPixelsDataFrame") %>% 
    as.data.frame() %>% 
    mutate(vims_dens = layer) %>% 
    dplyr::select(-layer)
  
  # predict from barrier model
  bar_mod_fc <- sdmTMB:::predict.sdmTMB(bar_mod_dens, newdata = bar_rast_df,
                                        return_tmb_object = TRUE)
  (bar_plt <- 
    plot_map(dat = bar_mod_fc$data, column = "exp(est)") +
    labs(fill = "Density\n(shoots m<sup>-2</sup>)",
         title = paste("Barrier -",yrs[yr])) +
    theme(legend.title = element_markdown(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "bottom")) 
  if (yr != 9){
    bar_index <- sdmTMB::get_index(bar_mod_fc, bias_correct = TRUE) %>%
      mutate(AIC = AIC(bar_mod_dens))
  } else {
    bar_index <- tibble(`_sdmTMB_time` = NA,
                         est = NA,
                         lwr = NA,
                         upr = NA,
                         log_est = NA,
                         se = NA,
                         max_gradient = NA,
                         bad_eig = NA,
                         AIC = NA)
  }
  
  # predict from stationary model
  mod_fc <- sdmTMB:::predict.sdmTMB(mod_dens, newdata = bar_rast_df,
                                        return_tmb_object = TRUE)
  (stat_plt <- 
    plot_map(dat = mod_fc$data, column = "exp(est)") +
    labs(fill = "Density\n(shoots m<sup>-2</sup>)",
         title = paste("Stationary -",yrs[yr])) +
    theme(legend.title = element_markdown(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "bottom"))
  if (yr != 9){
    stat_index <- sdmTMB::get_index(mod_fc, bias_correct = TRUE) %>%
      mutate(AIC = AIC(mod_dens))
  } else {
    stat_index <- tibble(`_sdmTMB_time` = NA,
                          est = NA,
                          lwr = NA,
                          upr = NA,
                          log_est = NA,
                          se = NA,
                          max_gradient = NA,
                          bad_eig = NA,
                          AIC = NA)
  }
  
  
  # predict from GAM
  mod_fc_gam <- 
    cbind(bar_rast_df, 
          predict.gam(mod_dens_gam, newdata = bar_rast_df,
                se.fit = T,
                type = "response")) %>% 
    dplyr::rename(est = fit)
  
  (gam_plt <- 
    plot_map(dat = mod_fc_gam, column = "est", sdmtmb = F) +
    labs(fill = "Density\n(shoots m<sup>-2</sup>)",
         title = paste("GAM -",yrs[yr])) +
    theme(legend.title = element_markdown(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "bottom"))
  
  assign("pred_shoot_dens", rbind(pred_shoot_dens,
                                  bind_rows(mod_fc_gam %>% 
                                              dplyr::select(x,y, est) %>% 
                                              as_tibble() %>% 
                                              mutate(year = yrs[yr],
                                                     model = "GAM"),
                                            mod_fc$data %>% 
                                              dplyr::select(x,y, est) %>% 
                                              mutate(year = yrs[yr],
                                                    model = "TMB-Stationary",
                                                    est = exp(est)),
                                            bar_mod_fc$data %>% 
                                              dplyr::select(x,y, est) %>% 
                                              mutate(year = yrs[yr],
                                                     model = "TMB-Barrier",
                                                     est = exp(est)))))
  
  assign("shoot_index", rbind(shoot_index,
                              bind_rows(stat_index %>% mutate(year = yrs[yr],
                                                              model = "TMB-Stationary"),
                                        bar_index %>% mutate(year = yrs[yr],
                                                            model = "TMB-Barrier"))))
  
  save(pred_shoot_dens,shoot_index, file = here::here("data/iid_mod_output.rdata"))
}

ggplot() +
  geom_raster(data = pred_shoot_dens %>% 
                filter(model != "TMB-Stationary") %>% 
                mutate(model = ifelse(model == "GAM", "Stationary","Barrier")),
                aes(x = x, y = y, fill = est)) +
  coord_sf()  +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom") +
  scale_fill_gradientn(colors = pals::ocean.deep(10)) +
  labs(fill = "Shoot density\n(shoots m<sup>-2</sup>)") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_markdown(),
        legend.position = "bottom") +
  facet_grid(model~year)


shoot_index %>% 
  add_row(year = rep(c(2013, 2014, 2016), each = 2),
          model = rep(c("TMB-Barrier",
                    "TMB-Stationary"), 3)) %>% 
ggplot() +
  geom_point(aes(y = est, x =year, color = model)) +
  geom_line(aes(y = est, x = year, color = model)) +
  geom_errorbar(aes(ymax = upr, ymin = lwr, x = year,
                  color = model),
              width = 0.1) +
  labs(y = "Shoots",
       title = "Relative eelgrass abundance",
       subtitle = "Hog Island Bay") + 
  theme(axis.title.x = element_blank()) +
  dream::theme_fade() +
  ggsci::scale_color_d3() +
  ggsci::scale_fill_d3()

