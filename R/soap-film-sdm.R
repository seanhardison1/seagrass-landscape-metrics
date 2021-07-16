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
library(dsm)
# Thanks to Gavin Simpson for providing the inspiration
# and a guide to fitting these models

# load data----
load(here::here("data/bound_sg.rdata"))
load(here::here("data/HogIslandBaySeagrass.rdata"))
load(here::here("data/processed_dem.rdata"))

hack <- 2500

# process----
train_poly1 <- hi_sg %>% 
  filter(YEAR %in% c(2017, 2018)) %>% 
  group_by(year = YEAR) %>% 
  dplyr::summarise() %>% 
  st_as_sf() %>% 
  smoothr::fill_holes(threshold = 50000) %>% 
  smoothr::drop_crumbs(threshold = 1e6)

tp_2017 <- train_poly1 %>% filter(year == 2017)
tp_2017$geometry <- tp_2017$geometry + hack
st_crs(tp_2017) <- st_crs(train_poly1)

train_poly <- 
  st_transform(tp_2017, st_crs(train_poly1)) %>% 
  bind_rows(.,train_poly1 %>% filter(year == 2018)) %>% 
  dplyr::select(-year)


ggplot() +
  geom_sf(data = train_poly)

syn_full1 <- 
  sg_bound %>% 
    filter(year %in% c(2017,2018), 
           str_detect(SiteName, "HI"),
           vims_pres == 1) %>%
    dplyr::select(shoots, year) %>% 
    mutate(year = factor(year, levels = c(2017, 2018)))

syn_full_2017 <- syn_full1 %>% filter(year == 2017)
syn_full_2017$geometry <- syn_full_2017$geometry + hack
st_crs(syn_full_2017) <- st_crs(syn_full1)

syn_full <- 
  st_transform(syn_full_2017, st_crs(syn_full1)) %>% 
  bind_rows(.,syn_full1 %>% filter(year == 2018)) %>% 
  arrange(desc(year)) %>% 
  # dplyr::select(-year) %>% 
  dream::sfc_as_cols() %>% 
  st_set_geometry(NULL)

syn_full1 %<>% 
  dream::sfc_as_cols() %>% 
  st_set_geometry(NULL)

syn_dat <- syn_full %>% slice(1:65)
head(syn_dat)
syn_dat_tprs <- syn_full1 %>% slice(1:65)
head(syn_dat_tprs)

tss_test_so  <- syn_full %>% slice(66:80)
tss_test_tprs  <- syn_full1 %>% slice(66:80)

ggplot() +
  geom_sf(data = train_poly) +
  geom_point(data = syn_dat, aes(x = x, y = y)) 

# geom_point(data = syn_dat, aes(x= x, y = y, color = shoots))

# sg_bound %>% 
#   filter(str_detect(SiteName, "HI")) %>% 
#   group_by(year, vims_pres) %>% 
#   dplyr::summarise(n = n())

# smooth with TPRSs----

# Use TPRS to model how shoot density changes spatially
syn_dat_tprs %<>% mutate(year = factor(year))
m1 <- gam(shoots ~ s(x, y, by = year), data = syn_dat_tprs, method = "REML")
summary(m1)

r <- raster(extent(hi_sg),
            res = 50,
            crs ="+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")

pdf1 <- fasterize(train_poly1 %>% filter(year == 2017), r) %>% 
  as("SpatialPixels") %>% 
  as.data.frame() %>% 
  mutate(year = 2017)

pdf2 <- fasterize(train_poly1 %>% filter(year == 2018), r) %>% 
  as("SpatialPixels") %>% 
  as.data.frame() %>% 
  mutate(year = 2018)

pdf <- bind_rows(pdf1, pdf2) %>% 
  mutate(year = factor(year))


##predictions
pdata <- pdf %>% mutate(shoots = predict(m1, newdata = pdf, type = "response"))
tmp <- pdata %>% 
  st_as_sf(coords = c("x","y"), 
           crs = st_crs(train_poly1)) %>% 
  st_intersection(train_poly1) %>% 
  dream::sfc_as_cols() %>% 
  st_set_geometry(NULL)

# visualize
ggplot() +
  geom_raster(data = tmp, aes(x = x, y = y, fill = shoots)) +
  geom_sf(data = train_poly1, fill = "transparent") +
  viridis::scale_fill_viridis(na.value = NA) +
  facet_wrap(~year) +
  geom_point(data = syn_dat_tprs,
             aes(x = x, y = y, size = shoots))

# soap-film smoothing----

# generate boundary for soap film smoother
crds <- st_coordinates(train_poly)[,c("X","Y","L2")]
crds1 <- crds[crds[, "L2"] == 1,]
crds2 <- crds[crds[, "L2"] == 2,]

crds <- 
  tibble(x = c(crds1[,1],crds2[,1]),
       y = c(crds1[,2],crds2[,2]),
       year = c(rep(2017, nrow(crds1)),
                rep(2018, nrow(crds2))))

names(crds) <- c("x", "y", "year")
bound <- split(crds, crds$year)
bound <-lapply(bound,`[`, c(1, 2))
nr <- seq(1, 2)
bound <-lapply(nr, function(n) as.list.data.frame(bound[[n]]))

## bound is a list of 3 lists (because 1 large polygon and 2 loops), each including 3 elements (utmx, utmy and f)
## an element f is set to a vector of zeros to set the boundary condition
bound <- lapply(nr, function(n) bound[[n]] <- c(bound[[n]], 
                                                list(f = rep(0, 
                                                             length(bound[[n]]$x)))))

# generate knots within the boundary
grid_out <- NULL
for (i in 1:nrow(train_poly)){
  grid_int <- 
    train_poly %>% 
    slice(i) %>% 
    st_make_grid(n=c(10,10), 
                 what = "centers", 
                 square=TRUE) %>% 
    st_intersection(train_poly)
  d = st_distance(grid_int, st_cast(train_poly,"MULTILINESTRING"))
  gridin2 = grid_int[d[,i]>units::set_units(0.075,km)]
  # plot(gridin2)
  assign('grid_out', rbind(grid_out, gridin2 %>% st_as_sf))
}


# plot(grid_out)

grid <- 
  st_coordinates(grid_out) %>% 
  as_tibble() %>% 
  dplyr::rename(x = X, y = Y)

ggplot() +
  geom_sf(data = train_poly) +
  # geom_point(data = syn_dat, aes(x = x, y = y)) +
  geom_point(data = grid,
             aes(x = x, y = y), alpha = 0.25)


# fit the model
m2 <- gam(shoots ~ s(x, y, 
                     bs = "so", 
                     xt = list(bnd = bound)),
          data = syn_dat, method = "REML", 
          knots = grid,
          na.action = na.omit)

summary(m2)

lims <- apply(crds, 2, range)
ylim <- lims[,2]
xlim <- lims[,1]

plot(m2, asp = 1, ylim = ylim, 
     xlim = xlim, se = FALSE, scheme = 2, main = "")

crds <- st_coordinates(train_poly)[,c("X","Y")]
grid.x <- with(m1$var.summary,
               seq(min(c(x, crds[,1])), max(c(x, crds[,1])), by = 50))
grid.y <- with(m1$var.summary,
               seq(min(c(y, crds[,2])), max(c(y, crds[,2])), by = 50))
pdata_grid <- with(m2$var.summary, expand.grid(x = grid.x, y = grid.y))
names(pdata_grid) <- c("x","y")

pdata2 <- transform(pdata_grid[, 1:2], shoots = predict(m2, newdata = pdata_grid,
                                                        type = "response"))

(soap <- 
  ggplot() +
  geom_raster(data = pdata2, aes(x = x, y = y, fill = shoots)) +
  geom_sf(data = train_poly, fill = "transparent") +
    geom_contour(data = pdata2, aes(x = x, y = y, z = shoots),
                 show.legend = F, color = "grey70",
                 bins = 6) +
    geom_text_contour(data = pdata2, aes(x = x, y = y, z = shoots),
                      size = 2, label.placement = label_placement_n(1),
                      skip = 1, color = "grey10") +
  viridis::scale_fill_viridis(na.value = NA,
                              limits = c(0,550)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 5),
        legend.title = element_textbox()) +
  labs(fill = "Shoots m<sup>-2</sup>"))

(tprs <- 
  ggplot() +
  geom_raster(data = tmp, aes(x = x, y = y, fill = shoots)) +
  geom_sf(data = train_poly1, fill = "transparent") +
  geom_contour(data = tmp, aes(x = x, y = y, z = shoots),
               show.legend = F, color = "grey70",
               bins = 8) +
  geom_text_contour(data = tmp, aes(x = x, y = y, z = shoots),
                          size = 2, label.placement = label_placement_n(1),
                    skip = 1, color = "grey10") +
  facet_wrap(~year) +
  viridis::scale_fill_viridis(na.value = NA,
                              limits = c(-250,600)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 5),
        legend.title = element_textbox()) +
  labs(subtitle = "TPRS smooth",
       fill = "Shoots m<sup>-2</sup>"))

soap + tprs + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave(filename = here::here("sdm_comparison.png"),
       dpi = 200, width = 7, height = 3)

tprs_test_pred <- predict(m1,newdata=tss_test_tprs, type = "response")
soap_test_pred <- predict(m2,newdata=tss_test_so, type = "response")

(oos_pred <- 
  AIC(m1, m2) %>% 
  mutate(RMSE = c(sqrt(mean((tss_test_tprs$shoots - tprs_test_pred)^2)),
                  sqrt(mean((tss_test_so$shoots - soap_test_pred)^2)))))

row.names(oos_pred) <- c("Kriging","Soap-film")

syn_full %>%
  mutate(soap_test_pred = as.numeric(predict(m2, newdata = ., type = "response")),
         tprs_test_pred = as.numeric(predict(m1, newdata = ., type = "response"))) %>% 
  gather(var, value, -x,-y,-shoots,-year) %>% 
  ggplot() +
    geom_point(aes(x = shoots, y = value)) +
    facet_wrap(.~var)
