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

# Get seagrass data 
load(here::here("data/bound_sg.rdata"))
load(here::here("data/vims_sg.rdata"))

coor_rs_km <- "+proj=utm +zone=18 +datum=NAD83 +units=km +no_defs"

sg_bound_figs <- sg_bound %>% 
  dplyr::rename(shoots = SHOOTS) %>% 
  mutate(meadow = case_when(meadow == "HI" ~ "Hog Island",
                            meadow == "SB" ~ "South Bay")) %>% 
  st_transform(crs = st_crs(coor_rs_km))

sg_bound_summ <- 
  sg_bound_figs %>% 
  group_by(year, meadow) %>% 
  dplyr::summarise(mean_shoot_dens = mean(shoots)) %>% 
  st_set_geometry(NULL)

vims_sg_figs <- vims_sg_proc %>% 
  filter(vims_dens > 0) %>% 
  mutate(year = factor(year, levels = 2007:2018)) %>% 
  st_transform(crs = st_crs(coor_rs_km))

hi <- st_intersection(vims_sg_figs,
                      st_read(here::here("data/HogIslandBay.kml")) %>% 
                        sf::st_transform(crs = st_crs(coor_rs_km)) %>% 
                        st_zm() %>% 
                        st_as_sf()) %>% 
  add_row(year = factor(c(2013, 2014, 2016)),
          vims_dens = c(NA, NA, NA)) 

sb <- st_intersection(vims_sg_figs,
              st_read(here::here("data/SouthBay.kml")) %>% 
                sf::st_transform(crs = st_crs(coor_rs_km)) %>% 
                st_zm() %>% 
                st_as_sf()) %>% 
  add_row(year = factor(c(2014, 2016)),
          vims_dens = c(NA, NA)) 

sg_area <- hi %>% 
  group_by(year) %>% 
  summarise() %>% 
  mutate(area = as.numeric(st_area(.)),
         meadow = "Hog Island") %>% 
  bind_rows(sb %>% 
              group_by(year) %>% 
              summarise() %>% 
              mutate(area = as.numeric(st_area(.)),
                     meadow = "South Bay")) %>% 
  mutate(area = ifelse(area == 0, NA, area)) %>% 
  left_join(.,sg_bound_summ)

ggplot(data = sg_area) +
  geom_point(aes(x = year, y = area, color = meadow))

ggplot(data = sg_area) +
  geom_point(aes(y = area, x = mean_shoot_dens)) +
  facet_wrap(~meadow)

(sg_dens_fig <- 
  ggplot(sg_bound_figs) +
  geom_boxplot(aes(x = factor(year), y = shoots, fill = meadow),
               position = position_dodge(),
               alpha = 0.5) +
  ggsci::scale_fill_d3() +
  dream::theme_fade() +
  labs(y = "Shoots m<sup>-2</sup>",
       fill = "Meadow") +
  theme(axis.title.y = element_markdown(),
        axis.title.x = element_blank()))

ggsave(sg_dens_fig, filename = here::here("figures/sg_dens_meadow_compare.png"),
       dpi = 300, width = 8, height = 5)


  ggplot() +
    geom_sf(data = vims_sg_figs) +
    facet_wrap(~year)
