library(tidyverse)
library(rvest)
library(sf)
library(rgdal)

site <-xml2::read_html("http://web.vims.edu/bio/sav/gis_data.html")
base <- "http://web.vims.edu/bio/sav/"


urls <- site %>% 
  html_nodes(., "a") %>% 
  html_attr(., "href") %>% 
  as_tibble() %>% 
  dplyr::filter(str_detect(value, "zip")) %>% 
  mutate(value = paste0(base, value)) %>% 
  pull(value)

for (i in 1:length(urls)){
  
  yr <- paste0("seagrass",str_extract(urls[i], "\\d+"),".zip")
  print(yr)
  download.file(url = urls[i], destfile = paste0(here::here("data/"),yr))
}

