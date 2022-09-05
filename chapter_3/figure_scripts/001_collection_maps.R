#!/usr/bin/env Rscript
# collection sites maps
# cjfiscus
# 2022-03-04

library(pacman)
p_load(readxl, ggplot2, ggmap, maps, mapdata,viridis, kgc, 
       cowplot, magrittr, maptools, raster, rgeos, rasterVis, viridis)

setwd("~/Desktop/FIGURE_DRAFTS_COWPEA/SCRIPTS/")

# import data
df<-read_excel("../data/IITA-Cowpea collection.xls")
df<-df[,c("Accession name", "Latitude", "Longitude", "Country of origin", "Biological status of accession", "In Core")]
df$`Accession name`<-toupper(df$`Accession name`)

samps<-read.table("../data/samples.txt")
samps$V1<-toupper(samps$V1)

## subset to samps
df<-df[df$`Accession name` %in% samps$V1,]

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
  geom_point(data=df, aes(x=Longitude, y=Latitude), color="blue", size=1, alpha=0.5) +
  theme_void() +
  coord_fixed(ratio=1)
ggsave("001A_collection_sites_world.pdf", p1, height=5, width=12)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3)+
  geom_point(data=df, aes(x=Longitude, y=Latitude), color="blue", size=1, alpha=0.5) +
  theme_void() +
  coord_fixed(xlim=c(-30, 60), ylim=c(-40, 40), ratio=1)
ggsave("001B_collection_sites_africa.pdf", p1, height=5, width=5)




p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3)+
  geom_point(data=df, aes(x=Longitude, y=Latitude, 
                          color=`Biological status of accession`), size=1, alpha=0.5) +
  theme_void() +
  coord_fixed(xlim=c(-30, 60), ylim=c(-40, 40), ratio=1) +
  theme(legend.position="top") +
  scale_color_brewer(palette="Set1")
ggsave("001C_collection_sites_africa_biostatus.pdf", p1, height=5, width=5)
