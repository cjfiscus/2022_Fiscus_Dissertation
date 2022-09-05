#!/usr/bin/env Rscript
# pca maps 
# cjfiscus
# 2022-03-07

library(pacman)
p_load(readxl, ggplot2, ggmap, maps, mapdata,viridis, kgc, 
       cowplot, magrittr, maptools, raster, rgeos, rasterVis, viridis)

setwd("~/Desktop/FIGURE_DRAFTS_COWPEA/SCRIPTS/")

# import data
df<-read_excel("../data/IITA-Cowpea collection.xls")
df<-df[,c("Accession name", "Latitude", "Longitude", "Country of origin", "Biological status of accession")]
names(df)[1]<-"sample"

samps<-read.table("../data/samples.txt")

## subset to samps
df<-df[df$sample %in% samps$V1,]

df1<-read.delim("../data/scores_pc.txt")
df1$sample<-row.names(df1)

df<-merge(df, df1, by="sample")

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

## world maps
p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3)+
  geom_point(data=df, aes(x=Longitude, y=Latitude, color=Axis1), size=1, alpha=0.5) +
  theme_void() +
  coord_fixed(ratio=1) +
  scale_color_viridis(option="turbo") +
  theme(legend.position="top")
ggsave("013A_pca_pc1_world.pdf", p1, height=5, width=12)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3)+
  geom_point(data=df, aes(x=Longitude, y=Latitude, color=Axis2), size=1, alpha=0.5) +
  theme_void() +
  coord_fixed(ratio=1) +
  scale_color_viridis(option="turbo") +
  theme(legend.position="top")
ggsave("013C_pca_pc2_world.pdf", p1, height=5, width=12)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3)+
  geom_point(data=df, aes(x=Longitude, y=Latitude, color=Axis3), size=1, alpha=0.5) +
  theme_void() +
  coord_fixed(ratio=1) +
  scale_color_viridis(option="turbo") +
  theme(legend.position="top")
ggsave("013E_pca_pc3_world.pdf", p1, height=5, width=12)

## africa maps
p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3)+
  geom_point(data=df, aes(x=Longitude, y=Latitude, color=Axis1), size=1, alpha=0.5) +
  theme_void() +
  coord_fixed(xlim=c(-30, 60), ylim=c(-40, 40), ratio=1) +
  scale_color_viridis(option="turbo") +
  theme(legend.position="top")
ggsave("013B_pca_pc1_africa.pdf", p1, height=5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3)+
  geom_point(data=df, aes(x=Longitude, y=Latitude, color=Axis2), size=1, alpha=0.5) +
  theme_void() +
  coord_fixed(xlim=c(-30, 60), ylim=c(-40, 40), ratio=1) +
  scale_color_viridis(option="turbo") +
  theme(legend.position="top")
ggsave("013D_pca_pc2_africa.pdf", p1, height=5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3)+
  geom_point(data=df, aes(x=Longitude, y=Latitude, color=Axis3), size=1, alpha=0.5) +
  theme_void() +
  coord_fixed(xlim=c(-30, 60), ylim=c(-40, 40), ratio=1) +
  scale_color_viridis(option="turbo") +
  theme(legend.position="top")
ggsave("013F_pca_pc3_africa.pdf", p1, height=5, width=5)