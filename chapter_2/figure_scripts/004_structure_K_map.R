#!/usr/bin/env Rscript
# collection sites maps
# cjfiscus
# 2022-06-29

library(pacman)
p_load(readxl, ggplot2, ggmap, maps, mapdata,viridis, kgc, 
       cowplot, magrittr, maptools, raster, rgeos, rasterVis, viridis, 
       scatterpie, tidyr, ggsci, ggsci)

setwd("~/Desktop/MANUSCRIPTS/FIGURE_DRAFTS_CAPSELLA/SCRIPTS/")

# import data
df<-read.delim("../data/rad_data_samples_cw.txt")
names(df)[1]<-"sample"

## import cluster assignments
df1<-read.delim("../data/cluster_assignments.txt")

df1<-merge(df1, df, by="sample")

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

## k = 2
df2<-df1[df1$K==2,]
df2<-df2[,c("sample", "Long", "Lat", "variable", "value", "assignment")]
df2<-spread(df2, variable, value)
df2<-na.omit(df2)
df2$assignment<-factor(df2$assignment, levels=c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Admixed"))

pal<-c("#4DBBD5FF", "#E64B35FF", "#00A087FF", "#3C5488FF", "dimgrey")
names(pal)<-c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Admixed")

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
                  aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal) +
  geom_hline(aes(yintercept=42.35)) +
  geom_hline(aes(yintercept=34.65))
ggsave("004A_structure_map_world_k2_marked.pdf", p1, height=5, width=12)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(-20, 60), ylim=c(30, 70), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004B_structure_map_weurope_k2.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(30, 160), ylim=c(10, 70), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004C_structure_map_asia_k2.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(110, 160), ylim=c(-45, -10), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004D_structure_map_oceania_k2.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(-85, -30), ylim=c(-60, 15), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004E_structure_map_samerica_k2.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(-180, -30), ylim=c(10, 80), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004F_structure_map_namerica_k2.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(-125, -110), ylim=c(30, 43), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004G_structure_map_cali_k2.pdf", p1, height=5.5, width=5)

### plot dist of clusters by latitude
quantile(test$Lat)

ggplot(df2, aes(x=assignment, y=Lat)) + geom_violin() +
  geom_hline(aes(yintercept=42))

##########
## k = 3
df2<-df1[df1$K==3,]
df2<-df2[,c("sample", "Long", "Lat", "variable", "value", "assignment")]
df2<-spread(df2, variable, value)
df2<-na.omit(df2)
df2$assignment<-factor(df2$assignment, levels=c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Admixed"))

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004A_structure_map_world_k3.pdf", p1, height=5, width=12)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(-20, 60), ylim=c(30, 70), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004B_structure_map_weurope_k3.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(30, 160), ylim=c(10, 70), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004C_structure_map_asia_k3.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(110, 160), ylim=c(-45, -10), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004D_structure_map_oceania_k3.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(-85, -30), ylim=c(-60, 15), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004E_structure_map_samerica_k3.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(-180, -30), ylim=c(10, 80), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004F_structure_map_namerica_k3.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(-125, -110), ylim=c(30, 43), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004G_structure_map_cali_k3.pdf", p1, height=5.5, width=5)

## k = 4
df2<-df1[df1$K==4,]
df2<-df2[,c("sample", "Long", "Lat", "variable", "value", "assignment")]
df2<-spread(df2, variable, value)
df2<-na.omit(df2)
df2$assignment<-factor(df2$assignment, levels=c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Admixed"))

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004A_structure_map_world_k4.pdf", p1, height=5, width=12)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(-20, 60), ylim=c(30, 70), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004B_structure_map_weurope_k4.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(30, 160), ylim=c(10, 70), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004C_structure_map_asia_k4.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(110, 160), ylim=c(-45, -10), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004D_structure_map_oceania_k4.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(-85, -30), ylim=c(-60, 15), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004E_structure_map_samerica_k4.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(-180, -30), ylim=c(10, 80), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004F_structure_map_namerica_k4.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(xlim=c(-125, -110), ylim=c(30, 43), ratio=1) +
  theme_void() +
  theme(legend.position="top") +
  scale_color_manual(values=pal)
ggsave("004G_structure_map_cali_k4.pdf", p1, height=5.5, width=5)
