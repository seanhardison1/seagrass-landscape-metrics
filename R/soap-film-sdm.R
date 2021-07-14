library(tidyverse)
library(sf)
library(raster)
library(sdmTMB)
library(rsample) 
library(patchwork)
library(INLA)
library(ggtext)
library(fasterize)
library(metR)
library(mgcv)
library(magrittr)

# Thanks to Gavin Simpson for providing the inspiration
# and a guide to fitting these models

# load data----
load(here::here("data/bound_sg.rdata"))
load(here::here("data/HogIslandBaySeagrass.rdata"))
load(here::here("data/processed_dem.rdata"))

# process----
train_poly <- hi_sg %>% 
  filter(YEAR == 2018) %>% 
  slice(1, 2, 5) %>% 
  st_union() %>% 
  st_as_sf() %>% 
  smoothr::fill_holes(threshold = 50000)

syn_full <- 
  sg_bound %>% 
    filter(year %in% 2018, str_detect(SiteName, "HI"),
           vims_pres == 1) %>%
    dplyr::select(shoots) %>% 
  dream::sfc_as_cols() %>% 
  st_set_geometry(NULL)

tss_split <- initial_split(syn_full, prop = .8)
syn_dat <- training(tss_split)
tss_test  <- testing(tss_split)

ggplot() +
  geom_sf(data = train_poly) +
  labs(title = "Hog Island meadow (2018)")
  # geom_point(data = syn_dat, aes(x= x, y = y, color = shoots))

# sg_bound %>% 
#   filter(str_detect(SiteName, "HI")) %>% 
#   group_by(year, vims_pres) %>% 
#   dplyr::summarise(n = n())

# smooth with TPRSs----

# Use TPRS to model how shoot density changes spatially
crds <- st_coordinates(train_poly)[,c("X","Y")]
m1 <- gam(shoots ~ s(x, y), data = syn_dat, method = "REML")
summary(m1)

grid.x <- with(m1$var.summary,
               seq(min(c(x, crds[,1])), max(c(x, crds[,1])), by = 50))
grid.y <- with(m1$var.summary,
               seq(min(c(y, crds[,2])), max(c(y, crds[,2])), by = 50))
pdata <- with(m1$var.summary, expand.grid(x = grid.x, y = grid.y))
names(pdata) <- c("x","y")

##predictions
pdata <- transform(pdata, shoots = predict(m1, pdata, type = "response"))
tmp <- pdata %>% 
  st_as_sf(coords = c("x","y"), 
           crs = st_crs(train_poly)) %>% 
  st_intersection(train_poly) %>% 
  dream::sfc_as_cols() %>% 
  st_set_geometry(NULL)

# visualize
ggplot() +
  geom_raster(data = tmp, aes(x = x, y = y, fill = shoots)) +
  geom_sf(data = train_poly, fill = "transparent") +
  viridis::scale_fill_viridis(na.value = NA)
  
# soap-film smoothing----

# generate boundary for soap film smoother
bound <- list(list(x = crds[,1], y = crds[,2], f = rep(0, nrow(crds))))

# generate knots within the boundary
grid <- st_make_grid(train_poly, n=c(8,8), 
                     what = "centers", 
                     square=TRUE) %>% 
  st_intersection(train_poly)  

d <- st_distance(grid, st_cast(train_poly,"MULTILINESTRING"))

# make sure knots aren't too close to the boundary
knots_fixed <- grid[d[,1]>units::set_units(0.025,km)] %>% 
  st_as_sf() %>% 
  dream::sfc_as_cols() %>% 
  st_as_sf() %>% 
  st_set_geometry(NULL) %>% 
  dplyr::rename(x = 1) %>% 
  as.matrix()

in_out <- in.out(bound, knots_fixed)

knots_fixed2 <- knots_fixed %>% 
  as_tibble() %>% 
  filter(in_out) %>% 
  as_tibble() #%>% 
  # dplyr::add_row(x = c(435348.9, 435204.5,435998.4),
  #          y = c(4140957, 4141215,4141370))

ggplot() +
  geom_sf(data = train_poly) +
  geom_point(data = knots_fixed2, 
             aes(x = x, y = y))

# fit the model
m2 <- gam(shoots ~ s(x, y, bs = "so", 
                     xt = list(bnd = bound),
                     k = 5),
          data = syn_dat, method = "REML", knots = knots_fixed2)

summary(m2)

lims <- apply(crds, 2, range)
ylim <- lims[,2]
xlim <- lims[,1]

plot(m2, asp = 1, ylim = ylim, 
     xlim = xlim, se = FALSE, scheme = 2, main = "")

pdata2 <- transform(pdata[, 1:2], shoots = predict(m2, newdata = pdata))

(soap <- 
  ggplot() +
  geom_raster(data = pdata2, aes(x = x, y = y, fill = shoots)) +
  geom_sf(data = train_poly, fill = "transparent") +
    geom_contour(data = pdata2, aes(x = x, y = y, z = shoots),
                 show.legend = F, color = "grey70",
                 bins = 8) +
    geom_text_contour(data = pdata2, aes(x = x, y = y, z = shoots),
                      size = 2, label.placement = label_placement_n(1),
                      skip = 1, color = "grey10") +
  viridis::scale_fill_viridis(na.value = NA,
                              limits = c(0,500)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 5),
        legend.title = element_textbox()) +
  labs(subtitle = "Soap-film smooth",
       title = "Comparing eelgrass SDM methods",
       fill = "Shoots m<sup>-2</sup>"))

(tprs <- 
  ggplot() +
  geom_raster(data = tmp, aes(x = x, y = y, fill = shoots)) +
  geom_sf(data = train_poly, fill = "transparent") +
  geom_contour(data = tmp, aes(x = x, y = y, z = shoots),
               show.legend = F, color = "grey70",
               bins = 8) +
  geom_text_contour(data = pdata2, aes(x = x, y = y, z = shoots),
                          size = 2, label.placement = label_placement_n(1),
                    skip = 1, color = "grey10") +
  viridis::scale_fill_viridis(na.value = NA,
                              limits = c(0,500)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 5),
        legend.title = element_textbox()) +
  labs(subtitle = "TPRS smooth",
       fill = "Shoots m<sup>-2</sup>"))

soap + tprs + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave(filename = here::here("sdm_comparison.png"),
       dpi = 200)

tprs_test_pred <- predict(m1,newdata=tss_test, type = "response")
soap_test_pred <- predict(m2,newdata=tss_test, type = "response")

(oos_pred <- 
  AIC(m1, m2) %>% 
  mutate(RMSE = c(sqrt(mean((tss_test$shoots - tprs_test_pred)^2)),
                  sqrt(mean((tss_test$shoots - soap_test_pred)^2)))))

row.names(oos_pred) <- c("Kriging","Soap-film")

syn_full %>%
  mutate(tprs_test_pred = predict(m1, newdata = .),
         soap_test_pred = predict(m2, newdata = .)) %>% 
  gather(var, value, -x,-y,-shoots) %>% 
  ggplot() +
    geom_point(aes(x = shoots, y = value)) +
    facet_wrap(.~var)
