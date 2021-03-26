library(tidyverse)
library(sf)
library(magrittr)
sg_bmass_canhei <- read_csv(here::here("data/sg_biomass_and_canopy_height.csv"),
                skip = 21) 
names(sg_bmass_canhei) <- str_to_lower(names(sg_bmass_canhei))
sg_bmass_canhei %<>% dplyr::rename(SiteName = plot)
save(sg_bmass_canhei, file = here::here("data/sg_bmass_canhei.rdata"))
