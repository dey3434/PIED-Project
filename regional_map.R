# map of PIED cored locations with species range distribution on it
# Author: Kelly Heilman

library(ggplot2)
library(sf)
library(tidyr)
library(tidyverse)
library(mapdata)
library(here)

# read in the PIED data used for modelling
PIED.all <- read.csv(here::here("data", "PIED_data.csv"))

PIED.ll <- unique(PIED.all[,c("LAT", "LON", "ELEV", "PLOT", "STATECD")])

# get the states using mapdata and maptools
# make a map of all of these:
all_states <- map_data("state")
states <- subset(all_states, region %in% c( "arizona", "utah", "new mexico", "colorado","idaho", "wyoming", "montana", "nevada", 
                                            "california", "oregon", "washington", "texas", "kansas", 
                                            "nebraska", "north dakota", "south dakota", "oklahoma") )
canada <- map_data("worldHires", "Canada")
mexico <- map_data("worldHires", "Mexico")


# note, I cloned this handy repository that has all of the Little species distribution maps
# to make this just update the dir.path to your cloned repo here: 
# https://github.com/wpetry/USTreeAtlas

dir.path = "/Users/kellyheilman/USTreeAtlas/shp/" # where the distribution shape files are

# the species in the little maps are the 1st 4 characters of the genus and species names
spp <- "pinuedul"
spp.distribution <- st_read(paste0(dir.path, spp, "/"))


if(file.exists(paste0(dir.path, spp, "/"))){
  spp.distribution <- st_read(paste0(dir.path, spp, "/"))
  
  #plt.spp <- plt.shp %>% filter(shp.code %in% spp)
  
  # plot distribution 
  spp.distribution %>% 
    ggplot() +
    geom_polygon(data = states, 
                 aes(x=long, y=lat, group = group), 
                 color = "black", fill = "white") +
    geom_polygon(data = mexico, 
                 aes(x=long, y=lat, group = group), 
                 color = "black", fill = "white") +
    geom_polygon(data = canada, 
                 aes(x=long, y=lat, group = group), 
                 color = "black", fill = "white") +
    geom_sf(alpha = 0.5, aes(fill = as.character(CODE)))+
    scale_fill_manual(values = c("1" = "forestgreen", "0" = "white"))+
    geom_point(data = PIED.ll, aes(x = LON, y = LAT), size = 0.75)+theme_bw()+
    coord_sf(xlim = c(-115, -103), ylim = c(30, 42))+theme(axis.title = element_blank(), legend.position = "none")
  
  ggsave(height = 6, width = 8, units = "in", here("images/", paste0(spp, "_distribution_map_zoom.png")))
  
  spp.distribution %>% 
    ggplot() +
    geom_polygon(data = states, 
                 aes(x=long, y=lat, group = group), 
                 color = "black", fill = "white") +
    geom_polygon(data = mexico, 
                 aes(x=long, y=lat, group = group), 
                 color = "black", fill = "white") +
    geom_polygon(data = canada, 
                 aes(x=long, y=lat, group = group), 
                 color = "black", fill = "white") +
    geom_sf(alpha = 0.5, aes(fill = as.character(CODE)))+
    scale_fill_manual(values = c("1" = "forestgreen", "0" = "white"))+
    geom_point(data = PIED.ll, aes(x = LON, y = LAT), size = 0.75)+theme_bw()+
    coord_sf(xlim = c(-120, -100), ylim = c(25, 48))+theme(axis.title = element_blank(), legend.position = "none")
  
  ggsave(height = 6, width = 8, units = "in", here("images/", paste0(spp, "_distribution_map_full.png")))
  
  
}  

