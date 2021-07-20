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


# helper function
sfc_as_cols <- function(x, geometry, names = c("x","y")) {
  if (missing(geometry)) {
    geometry <- sf::st_geometry(x)
  } else {
    geometry <- rlang::eval_tidy(enquo(geometry), x)
  }
  stopifnot(inherits(x,"sf") && inherits(geometry,"sfc_POINT"))
  ret <- sf::st_coordinates(geometry)
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[ , !names(x) %in% names]
  ret <- setNames(ret,names)
  dplyr::bind_cols(x,ret)
}

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
  sfc_as_cols() %>% 
  st_set_geometry(NULL) %>% 
  mutate(x = x/1000,
         y = y/1000)

syn_full1 %<>% 
  sfc_as_cols() %>% 
  st_set_geometry(NULL)%>% 
  mutate(x = x/1000,
         y = y/1000)

syn_dat <- syn_full %>% slice(1:80)
head(syn_dat)
syn_dat_tprs <- syn_full1 %>% slice(1:80)
head(syn_dat_tprs)

# 
# ggplot() +
#   geom_sf(data = train_poly) +
#   geom_point(data = syn_dat, aes(x = x, y = y)) 

# geom_point(data = syn_dat, aes(x= x, y = y, color = shoots))

# sg_bound %>% 
#   filter(str_detect(SiteName, "HI")) %>% 
#   group_by(year, vims_pres) %>% 
#   dplyr::summarise(n = n())

# smooth with TPRSs----

# Use TPRS to model how shoot density changes spatially
syn_dat_tprs %<>% mutate(year = factor(year))
m1 <- gam(shoots ~ s(x, y, by = year, bs = "gp"), data = syn_dat_tprs, method = "REML",
          family = tw(link = "log"))
summary(m1)
train_poly1 %<>% mutate(geometry = geometry/1000)

r <- raster(extent(train_poly1),
            res = 50/1000)

pdf1 <- fasterize(train_poly1 %>% filter(year == 2017), r) %>% 
  as("SpatialPixels") %>%
  as.data.frame() %>% 
  mutate(year = 2017)

pdf2 <- fasterize(train_poly1 %>% filter(year == 2018), r) %>% 
  as("SpatialPixels") %>% 
  as.data.frame() %>% 
  mutate(year = 2018)

pdf <- bind_rows(pdf1, pdf2) %>% 
  mutate(year = factor(year)) #%>% 
  # mutate(x = x/1000,
         # y = y/1000)

##predictions
pdata <- pdf %>% mutate(shoots = predict(m1, newdata = pdf, 
                                         type = "response"))
tmp <- pdata %>% 
  st_as_sf(coords = c("x","y"), 
           crs = st_crs(train_poly1)) %>% 
  st_intersection(train_poly1) %>% 
  sfc_as_cols() %>% 
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
                rep(2018, nrow(crds2)))) %>% 
  mutate(x = x/1000,
         y = y/1000)

names(crds) <- c("x", "y", "year")
bound <- split(crds, crds$year)
bound <-lapply(bound,`[`, c(1, 2))
nr <- seq(1, 2)
bound <-lapply(nr, function(n) as.list.data.frame(bound[[n]]))

# # set boundary of smooth to zero density - if not set, mgcv will estimate the boundary
# bound <- lapply(nr, function(n) bound[[n]] <- c(bound[[n]],
#                                                 list(f = rep(200,
#                                                              length(bound[[n]]$x)))))

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
  dplyr::rename(x = X, y = Y) %>% 
  mutate(x = x/1000,
         y = y/1000) %>% 
  add_row(x = c(435.3, 437.7, 435.2, 437.82),
          y = c(4141.1, 4143.864, 4141.3, 4143.6))


ggplot() +
  geom_sf(data = train_poly %>% mutate(geometry = geometry/1000)) +
  # geom_point(data = syn_dat, aes(x = x, y = y)) +
  geom_point(data = grid,
             aes(x = x, y = y), alpha = 0.25)


# fit the model
m2 <- gam(shoots ~ s(year, bs = "re") + 
                     s(x, y,
                     bs = "so",
                     xt = list(bnd = bound)),
          data = syn_dat, method = "REML",
          knots = grid)

m3 <- gam(shoots ~  s(x, y,
              bs = "so",
              xt = list(bnd = bound)),
          data = syn_dat, method = "REML",
          knots = grid)


summary(m3)
AIC(m1, m2, m3)
gratia::draw(m2)
gratia::appraise(m3)

lims <- apply(crds, 2, range)
ylim <- lims[,2]
xlim <- lims[,1]

plot(m3, asp = 1, ylim = ylim, 
     xlim = xlim, se = FALSE, scheme = 2, main = "")

crds <- st_coordinates(train_poly)[,c("X","Y")]
crds[,1] <- crds[,1]/1000
crds[,2] <- crds[,2]/1000

grid.x <- with(m3$var.summary,
               seq(min(c(x, crds[,1])), max(c(x, crds[,1])), by = 50/1000))
grid.y <- with(m3$var.summary,
               seq(min(c(y, crds[,2])), max(c(y, crds[,2])), by = 50/1000))
pdata_grid <- with(m3$var.summary, expand.grid(x = grid.x, y = grid.y,
                                               year = c(2017, 2018)))
names(pdata_grid) <- c("x","y","year")

pdata2 <- transform(pdata_grid, shoots = predict(m2, newdata = pdata_grid,
                                                        type = "response"))

(soap <- 
  ggplot() +
  geom_raster(data = pdata2, aes(x = x, y = y, fill = shoots)) +
  geom_sf(data = train_poly %>% mutate(geometry = geometry/1000), fill = "transparent") +
    geom_contour(data = pdata2, aes(x = x, y = y, z = shoots),
                 show.legend = F, color = "grey70",
                 bins = 6) +
    geom_text_contour(data = pdata2, aes(x = x, y = y, z = shoots),
                      size = 2, label.placement = label_placement_n(1),
                      skip = 1, color = "grey10") +
  viridis::scale_fill_viridis(na.value = NA,
                              limits = c(0,600)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 10),
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
                              limits = c(0,800)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 5),
        legend.title = element_textbox()) +
  labs(subtitle = "TPRS smooth",
       fill = "Shoots m<sup>-2</sup>"))

soap + tprs + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave(filename = here::here("sdm_comparison.png"),
       dpi = 200, width = 7, height = 3)

tprs_test_pred2 <- predict(m1,newdata=syn_dat_tprs, type = "response")
soap_test_pred2 <- predict(m3,newdata=syn_dat, type = "response")

(is_pred <- 
    AIC(m1, m3) %>% 
    mutate(RMSE = c(sqrt(mean((syn_dat_tprs$shoots - tprs_test_pred2)^2)),
                    sqrt(mean((syn_dat$shoots - soap_test_pred2)^2)))))

row.names(oos_pred) <- c("Kriging","Soap-film")

sf_pred <- syn_dat %>% 
  mutate(soap_test_pred = as.numeric(predict(m3, newdata = ., type = "response")))

tprs_pred <- syn_dat_tprs %>% 
  mutate(tprs_test_pred = as.numeric(predict(m1, newdata = ., type = "response")))

bind_rows(sf_pred, tprs_pred) %>%
  gather(var, value, -x,-y,-shoots,-year) %>% 
  ggplot() +
    geom_point(aes(x = shoots, y = value)) +
    facet_wrap(year~var, scales = "free") +
    geom_abline(aes(intercept = 0, slope = 1))
