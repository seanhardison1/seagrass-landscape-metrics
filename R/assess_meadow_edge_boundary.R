library(tidyverse)
library(sf)
library(raster)
library(caret)
library(patchwork)
library(INLA)
library(ggtext)
library(fasterize)

# Get seagrass data from presentation-------------------------
load(here::here("data/bound_sg.rdata"))

sg_bound %>% 
  mutate(vims_dens = factor(vims_dens)) %>% 
  filter(!is.na(vims_dens)) %>% 
  ggplot() +
    geom_vline(aes(xintercept = 0), color = "grey20") +
  geom_text(aes(x=70, label="\nOutside the meadow", y=550), 
            colour="brown", angle=-90, text=element_text(size=11), check_overlap = T) +
  geom_text(aes(x=-10, label="Inside the meadow\n", y=550), 
            colour="darkgreen", angle=90, text=element_text(size=11), check_overlap = T) +
    geom_point(aes(x = edge_dist,y = shoots, color = vims_dens)) +
    ggsci::scale_color_d3() +
    theme_minimal() +
    labs(y = "Known shoot density m<sup>-2</sup>",
         x = "Estimated distance to meadow edge (m)",
         color = "VIMS imagery\n density classes",
         title = "How well does the VIMS imagery capture the meadow edge?") +
    theme(legend.position = c(0.9,0.8),
          axis.title.y = element_markdown(size = 14),
          title = element_text(size = 14),
          axis.title.x = element_text(size = 14))
  

ggsave(filename = here::here("data/meadow_edge_assess.png"),device = "png", dpi = 300)
