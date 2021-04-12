library(tidyverse)
library(sf)
library(tsibble)
library(forecast)
library(fable)
load(here::here("data/seagrass_cover.rdata"))
load(here::here("data/seagrass_area.rdata"))

seagrass_coverage_interp <- 
  seagrass_coverage %>%
    tsibble(index = year, key = c(season, site)) %>%
    mutate(meadow = str_extract(site, "(?:(?!_).)*"),
           site = as.character(site)) %>% 
    fill_gaps() %>% 
    mutate(across(contains('sg_cover'), 
                  .fns = list(interp = ~na.interp(., linear = T))),
           interp = ifelse(is.na(sg_cover_1000), TRUE, FALSE))

ggplot(data = seagrass_coverage_interp) +
  geom_line(aes(x = year, y = sg_cover_1000_interp, group = site), color = "purple") +
  geom_point(aes(x = year, y = sg_cover_1000_interp, group = site), color = "purple") +
  geom_line(aes(x = year, y = sg_cover_1000, color = meadow, group = site)) +
  geom_point(aes(x = year, y = sg_cover_1000, color = meadow, group = site)) +
  facet_wrap(.~season)

save(seagrass_coverage, seagrass_coverage_interp,
     file = here::here("data/seagrass_cover.rdata"))
