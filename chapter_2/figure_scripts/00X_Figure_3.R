#!/usr/bin/env Rscript
# Figure 3
# cjfiscus
# 2022-08-01

library(pacman)
p_load(readxl, ggplot2, cowplot, ggsci, tidyr, stringr, kgc,
       ggmap, maps, mapdata, maptools, raster, rgeos, rasterVis, magrittr,
       dplyr, patchwork, viridis)

setwd("~/Desktop/MANUSCRIPTS/FIGURE_DRAFTS_CAPSELLA/SCRIPTS/")

# import coords
df1<-read.delim("~/Desktop/MANUSCRIPTS/FIGURE_DRAFTS_CAPSELLA/DATA/rad_data_samples_cw.txt")
names(df1)[1]<-"sample"
df1<-df1[,c("sample", "Lat", "Long")]

# import dist
df<-read.delim("../data/mean_dist_groups.txt")

## prep world map
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

##########
lst<-c("1118-03-00-00",
"1895-04-00-00",
"2184-10-00-00",
"2229-06-00-00")

# A 
df2<-df[df$comp=="NA_CA",]
m<-merge(df2, df1, by="sample")
m<-m[!m$sample %in% lst,]

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=m, 
             aes(x=Long, y=Lat, color=mean_dist)) +
  coord_fixed(ratio=1) +
  theme_void() +
  theme(legend.position="right", plot.margin = unit(c(0, 0, 0, 0), "null")) +
  scale_color_viridis(option="H", direction=-1)

# B
df2<-df[df$comp=="NA_other",]
m<-merge(df2, df1, by="sample")
m<-m[!m$sample %in% lst,]

p2<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=m, 
             aes(x=Long, y=Lat, color=mean_dist)) +
  coord_fixed(ratio=1) +
  theme_void() +
  theme(legend.position="right", plot.margin = unit(c(0, 0, 0, 0), "null")) +
  scale_color_viridis(option="H", direction=-1)

# C 
df2<-df[df$comp=="Japan",]
m<-merge(df2, df1, by="sample")
m<-m[!m$sample %in% lst,]

p3<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=m, 
             aes(x=Long, y=Lat, color=mean_dist)) +
  coord_fixed(ratio=1) +
  theme_void() +
  theme(legend.position="right", plot.margin = unit(c(0, 0, 0, 0), "null")) +
  scale_color_viridis(option="H", direction=-1)

# D 
df2<-df[df$comp=="GB_Ireland",]
m<-merge(df2, df1, by="sample")
m<-m[!m$sample %in% lst,]

p4<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=m, 
             aes(x=Long, y=Lat, color=mean_dist)) +
  coord_fixed(ratio=1) +
  theme_void() +
  theme(legend.position="right", plot.margin = unit(c(0, 0, 0, 0), "null")) +
  scale_color_viridis(option="H", direction=-1)
##########
# panels
ggsave("~/Desktop/A.pdf", p1, height=3, width=6)
ggsave("~/Desktop/B.pdf", p2, height=3, width=6)
ggsave("~/Desktop/C.pdf", p3, height=3, width=6)
ggsave("~/Desktop/D.pdf", p4, height=3, width=6)


leg<-get_legend(p1 + theme(legend.position="top"))
ggsave("~/Desktop/leg.pdf",leg, height=2, width=2)
