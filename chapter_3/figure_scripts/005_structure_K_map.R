#!/usr/bin/env Rscript
# collection sites maps
# cjfiscus
# 2022-03-04

library(pacman)
p_load(readxl, ggplot2, ggmap, maps, mapdata,viridis, kgc, 
       cowplot, magrittr, maptools, raster, rgeos, rasterVis, viridis, 
       scatterpie, tidyr, ggsci)

setwd("~/Desktop/MANUSCRIPTS/FIGURE_DRAFTS_COWPEA/SCRIPTS/")

# import data
df<-read_excel("../data/IITA-Cowpea collection.xls")
df<-df[,c("Accession name", "Latitude", "Longitude", "Country of origin", "Biological status of accession", "Collecting/acquisition source")]
df$`Accession name`<-toupper(df$`Accession name`)
names(df)[1]<-"sample"
df<-as.data.frame(df)

samps<-read.table("../data/samples.txt")
samps$V1<-toupper(samps$V1)

## subset to samps
df<-df[df$sample %in% samps$V1,]

## subset to lines we can't reliably
df3<-df[df$`Biological status of accession`=="Breeding / Research material" |
          df$`Collecting/acquisition source`=="Market or shop" |
          df$`Collecting/acquisition source`=="Institute, Experimental station, Research",]

## import cluster assignments
df1<-read.delim("../data/cluster_assignments.txt")
df1$sample<-toupper(df1$sample)

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
df2<-merge(df2, df, by="sample")
df2<-df2[,c("sample", "Longitude", "Latitude", "variable", "value")]
df2<-spread(df2, variable, value)
df2<-na.omit(df2)

# remove lines we can't reliably place
df2<-df2[!df2$sample %in% df3$sample,]

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_scatterpie(data=df2, 
                  aes(x=Longitude, y=Latitude, group=sample), 
                  cols=c("Cluster1", "Cluster2"), pie_scale=0.5) +
  coord_fixed(ratio=1) +
  theme_void() +
  scale_fill_npg() +
  theme(legend.position="top")
ggsave("005A_structure_map_world_k2.pdf", p1, height=5, width=12)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_scatterpie(data=df2, 
                  aes(x=Longitude, y=Latitude, group=sample), 
                  cols=c("Cluster1", "Cluster2"), pie_scale=0.3) +
  coord_fixed(xlim=c(-20, 50), ylim=c(-40, 40), ratio=1) +
  theme_void() +
  scale_fill_npg() +
  theme(legend.position="top")
ggsave("005B_structure_map_africa_k2.pdf", p1, height=5.5, width=5)

df3<-df1[df1$K==2,]
df3<-merge(df3, df, by="sample")
df3<-df3[,c("sample", "Longitude", "Latitude", "variable", "value", "assignment")]
pal<-c("#D55E00", "#009E73", "#56B4E9", "#CC79A7", "dimgrey")
names(pal)<-c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Admixed")
p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df3, 
                  aes(x=Longitude, y=Latitude, color=assignment), alpha=0.5) +
  coord_fixed(xlim=c(-20, 50), ylim=c(-40, 40), ratio=1) +
  theme_void() +
  scale_color_manual(values=pal) +
  theme(legend.position="none")
ggsave("~/Desktop/k2_africa.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_scatterpie(data=df2, 
                  aes(x=Longitude, y=Latitude, group=sample), 
                  cols=c("Cluster1", "Cluster2"), pie_scale=0.1) +
  coord_fixed(xlim=c(-20, 20), ylim=c(0, 20), ratio=1) +
  theme_void() +
  scale_fill_npg() +
  theme(legend.position="top")
ggsave("005C_structure_map_westafrica_k2.pdf", p1, height=5.5, width=5)

##########
## k = 3
df2<-df1[df1$K==3,]

df2<-merge(df2, df, by="sample")
df2<-df2[,c("sample", "Longitude", "Latitude", "variable", "value")]
df2<-spread(df2, variable, value)
df2<-na.omit(df2)

# remove lines we can't reliably place
df2<-df2[!df2$sample %in% df3$sample,]

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_scatterpie(data=df2, 
                  aes(x=Longitude, y=Latitude, group=sample), 
                  cols=c("Cluster1", "Cluster2", "Cluster3"), pie_scale=0.5) +
  coord_fixed(ratio=1) +
  theme_void() +
  scale_fill_npg() +
  theme(legend.position="top")
ggsave("005D_structure_map_world_k3.pdf", p1, height=5, width=12)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_scatterpie(data=df2, 
                  aes(x=Longitude, y=Latitude, group=sample), 
                  cols=c("Cluster1", "Cluster2", "Cluster3"), pie_scale=0.3) +
  coord_fixed(xlim=c(-20, 50), ylim=c(-40, 40), ratio=1) +
  theme_void() +
  scale_fill_npg() +
  theme(legend.position="top")
ggsave("005E_structure_map_africa_k3.pdf", p1, height=5.5, width=5)

df3<-df1[df1$K==3,]
df3<-merge(df3, df, by="sample")
df3<-df3[,c("sample", "Longitude", "Latitude", "variable", "value", "assignment")]
pal<-c("#D55E00", "#009E73", "#56B4E9", "#CC79A7", "dimgrey")
names(pal)<-c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Admixed")
p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df3, 
             aes(x=Longitude, y=Latitude, color=assignment), alpha=0.5) +
  coord_fixed(xlim=c(-20, 50), ylim=c(-40, 40), ratio=1) +
  theme_void() +
  scale_color_manual(values=pal) +
  theme(legend.position="none")
ggsave("~/Desktop/k3_africa.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_scatterpie(data=df2, 
                  aes(x=Longitude, y=Latitude, group=sample), 
                  cols=c("Cluster1", "Cluster2", "Cluster3"), pie_scale=0.1) +
  coord_fixed(xlim=c(-20, 20), ylim=c(0, 20), ratio=1) +
  theme_void() +
  scale_fill_npg() +
  theme(legend.position="top")
ggsave("005F_structure_map_westafrica_k3.pdf", p1, height=5.5, width=5)

## k = 4
df2<-df1[df1$K==4,]

df2<-merge(df2, df, by="sample")
df2<-df2[,c("sample", "Longitude", "Latitude", "variable", "value")]
df2<-spread(df2, variable, value)
df2<-na.omit(df2)

# remove lines we can't reliably place
df2<-df2[!df2$sample %in% df3$sample,]

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_scatterpie(data=df2, 
                  aes(x=Longitude, y=Latitude, group=sample), 
                  cols=c("Cluster1", "Cluster2", "Cluster3", "Cluster4"), pie_scale=0.5) +
  coord_fixed(ratio=1) +
  theme_void() +
  scale_fill_npg() +
  theme(legend.position="top")
ggsave("005G_structure_map_world_k4.pdf", p1, height=5, width=12)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_scatterpie(data=df2, 
                  aes(x=Longitude, y=Latitude, group=sample), 
                  cols=c("Cluster1", "Cluster2", "Cluster3", "Cluster4"), pie_scale=0.3) +
  coord_fixed(xlim=c(-20, 50), ylim=c(-40, 40), ratio=1) +
  theme_void() +
  scale_fill_npg() +
  theme(legend.position="top")
ggsave("005H_structure_map_africa_k4.pdf", p1, height=5.5, width=5)

df3<-df1[df1$K==4,]
df3<-merge(df3, df, by="sample")
df3<-df3[,c("sample", "Longitude", "Latitude", "variable", "value", "assignment")]
pal<-c("#D55E00", "#009E73", "#56B4E9", "#CC79A7", "dimgrey")
names(pal)<-c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Admixed")
p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df3, 
             aes(x=Longitude, y=Latitude, color=assignment), alpha=0.5) +
  coord_fixed(xlim=c(-20, 50), ylim=c(-40, 40), ratio=1) +
  theme_void() +
  scale_color_manual(values=pal) +
  theme(legend.position="top")
ggsave("~/Desktop/k4_africa.pdf", p1, height=5.5, width=5)

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_scatterpie(data=df2, 
                  aes(x=Longitude, y=Latitude, group=sample), 
                  cols=c("Cluster1", "Cluster2", "Cluster3", "Cluster4"), pie_scale=0.1) +
  coord_fixed(xlim=c(-20, 20), ylim=c(0, 20), ratio=1) +
  theme_void() +
  scale_fill_npg() +
  theme(legend.position="top")
ggsave("005I_structure_map_westafrica_k4.pdf", p1, height=5.5, width=5)

##########
# map of regions of Africa colored
## add regional data
region<-read.csv("~/Desktop/FIGURE_DRAFTS_COWPEA/DATA/UNSD_M49_2022.csv")
region<-region[region$Region.Name=="Africa",]
region$area<-ifelse(region$Sub.region.Name=="Northern Africa", "Northern Africa", region$Intermediate.Region.Name)
region<-region[,c("area", "Country.or.Area")]
names(region)<-c("area", "region")
add<-as.data.frame(rbind(c("Western Africa", "Ivory Coast"),
                         c("Middle Africa", "Republic of Congo"),
                         c("Eastern Africa", "Tanzania")))

names(add)<-c("area", "region")
region<-as.data.frame(rbind(region, add))

world <- map_data("world")
world2<-world
world2<-merge(world2, region, by="region", all.x=T)
ggplot() +
  geom_map(
    data = world2, map = world,
    aes(long, lat, map_id = region,
        fill = area),
    color = "black", size = 0.1
  ) +
  coord_fixed(xlim=c(-20, 50), ylim=c(-40, 40), ratio=1) + theme_void()

