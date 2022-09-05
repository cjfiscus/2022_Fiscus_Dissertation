#!/usr/bin/env Rscript
# collection sites maps
# cjfiscus
# 2022-06-29

library(pacman)
p_load(readxl, ggplot2, ggmap, maps, mapdata,viridis, kgc, 
       cowplot, magrittr, maptools, raster, rgeos, rasterVis, viridis)

setwd("~/Desktop/MANUSCRIPTS/FIGURE_DRAFTS_CAPSELLA//SCRIPTS/")

# import data
df<-read.delim("../data/rad_data_samples_cw.txt")

### plot map
## produce map
map<-map_data("world")

#Defining a general CRS
mycrs <- "+proj=longlat +datum=WGS84 +no_defs"

#Using the original maps package, then converting map into SpatialPolygons object
world <- maps::map("world", fill=TRUE) %$% 
  maptools::map2SpatialPolygons(., IDs=names,proj4string=CRS(mycrs))

#The resulting map has self intersection problems so any further operation reports errors; using buffers of width 0 is a fast fix
while(rgeos::gIsValid(world)==FALSE){
  world <- rgeos::gBuffer(world, byid = TRUE, width = 0, quadsegs = 5, capStyle = "ROUND")
}

#Dissolving polygon's limits
world <- raster::aggregate(world)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3)+
  geom_point(data=df, aes(x=Long, y=Lat), color="blue", size=1, alpha=0.5) +
  theme_void() +
  coord_fixed(ratio=1)
ggsave("001A_collection_sites_world.pdf", p1, height=5, width=12)
