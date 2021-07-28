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
library(sp)
library(gstat)

# load data----
load(here::here("data/bound_sg.rdata"))
load(here::here("data/HogIslandBaySeagrass.rdata"))
load(here::here("data/processed_dem.rdata"))


# process----
y <- 2018
# seagrass meadow boundary
train_poly1 <- hi_sg %>% 
  filter(YEAR %in% y) %>% 
  group_by(year = YEAR) %>% 
  dplyr::summarise() %>% 
  st_as_sf() %>% 
  smoothr::fill_holes(threshold = 50000) %>% 
  smoothr::drop_crumbs(threshold = 1e6) 

# data as sf
syn_full1 <- 
  sg_bound %>% 
  filter(year %in% y, 
         str_detect(SiteName, "HI"),
         vims_pres == 1) %>%
  st_transform(., crs = st_crs(train_poly1)) %>% 
  mutate(sg_dist = as.numeric(st_distance(., 
                                          st_cast(train_poly1,
                                                  "MULTILINESTRING")))) %>% 
  dplyr::select(shoots, year, sg_dist) %>% 
  mutate(year = factor(year, levels = y))

# generate prediction grid----
r <- raster(extent(train_poly1),
            res = 50)

pdf1 <- fasterize(train_poly1, r) %>% 
  as("SpatialPixels") %>%
  as.data.frame() 

# prediction grid for kriging model----
pdf2 <- pdf1
gridded(pdf2) = ~x+y
crs(pdf2) <- "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs"

# generate boundary for soap film----
crds <- st_coordinates(train_poly1)[,c("X","Y","L2")]
crds <- 
  tibble(x = crds[,1],
         y = crds[,2])
bound <- list(as.list.data.frame(crds))

# knots----
grid_int <- 
  train_poly1 %>% 
  st_make_grid(n=c(10,10), 
               what = "centers", 
               square=TRUE) %>% 
  st_intersection(.,train_poly1)
d = st_distance(grid_int, st_cast(train_poly1,"MULTILINESTRING"))
grid_out = grid_int[d>units::set_units(0.075,km)]
grid <- 
  st_coordinates(grid_out) %>% 
  as_tibble() %>% 
  dplyr::rename(x = X, y = Y) %>% 
  add_row(x = c(435300,  435200),
          y = c(4141100,  4141300))

# ggplot(grid) +
#   geom_point(aes(x = x, y = y)) +
#   geom_sf(data = train_poly1, fill = "transparent")

# plot(bound)
# points(grid$x, grid$y)
# full data as df----
syn_full_df <- syn_full1 %>%   
  dream::sfc_as_cols() %>% 
  st_set_geometry(NULL)

# 10-fold CV----
yourData<-syn_full_df[sample(nrow(syn_full_df)),]

#Create 10 equally size folds
folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)

# fit models and get summary stats
set.seed(2012)
pred_summs <- NULL
for(i in 1:10){
  message(i)
  #data segmentation
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- yourData[testIndexes, ]
  trainData <- yourData[-testIndexes, ]
  
  # soap film
  soap_mod <- gam(shoots ~  s(x, y,
                              bs = "so", k = 7,
                              xt = list(bnd = bound, nmax = 500)),
                  data = trainData, method = "REML",
                  knots = grid)
  # gratia::draw(soap_mod)
  # gratia::appraise(soap_mod)
  soap_pred <- predict(soap_mod, newdata = testData, type = "response")
  soap_rmse <- sqrt(mean((soap_pred - testData$shoots)^2))
  
  # TPRS
  tprs_mod <- gam(shoots ~  s(x, y),
                  data = trainData, method = "REML")
  tprs_pred <- predict(tprs_mod, newdata = testData)
  tprs_rmse <- sqrt(mean((tprs_pred - testData$shoots)^2))
  
  # data for kriging model----
  # fit variogram
  kg_df <- trainData
  coordinates(kg_df) = ~x+y
  crs(kg_df) <- "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs"
  
  kg_test <- testData
  coordinates(kg_test) <- ~x+y
  crs(kg_test) <- "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs"
  
  # fit variogram (ordinary krigin)
  v <- variogram(shoots~1, data = kg_df)
  vfit <- fit.variogram(v, vgm("Gau"))
  # print(plot(v, vfit))
  
  # ordinary kriging----
  ok_mod <- krige(formula = shoots~1, 
             kg_df, 
             newdata = kg_test, 
             model = vfit)
  # spplot(x["var1.pred"], main = "ordinary kriging predictions")
  ok_pred_test <- ok_mod@data$var1.pred
  ok_rmse <- sqrt(mean((ok_pred_test - testData$shoots)^2))
  
  # fit variogram (universal krigin)
  v2 <- variogram(shoots~sg_dist, data = kg_df)
  vfit2 <- fit.variogram(v2, vgm("Mat"), fit.method = 6)
  print(plot(v2, vfit2))
  
  # universal kriging----
  uni_mod <- krige(formula = shoots~sg_dist, 
                  kg_df, 
                  newdata = kg_test, 
                  model = vfit2)
  # spplot(x["var1.pred"], main = "ordinary kriging predictions")
  uni_pred_test <- uni_mod@data$var1.pred
  uni_rmse <- sqrt(mean((uni_pred_test - testData$shoots)^2))
  
  assign("pred_summs",rbind(
    tibble(ok_rmse,
         soap_rmse,
         tprs_rmse,
         uni_rmse,
         fold = i),
    pred_summs))
}
# visualize outputs
pred_summs %>% 
  gather(var, value,-fold) %>% 
ggplot() +
  geom_violin(aes(x = var, y = value),
              draw_quantiles = c(0.25, 0.5, 0.75))

pred_summs %>% 
  gather(var, value,-fold) %>% 
  group_by(var) %>% 
  dplyr::summarise(pred_rmse = mean(value))
