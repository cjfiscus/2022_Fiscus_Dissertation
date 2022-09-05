#!/usr/bin/env Rscript
# structure barplots
# cjfiscus
# 2022-06-30

library(pacman)
p_load(readxl, ggplot2, cowplot, ggsci, tidyr, stringr, kgc)

setwd("~/Desktop/MANUSCRIPTS/FIGURE_DRAFTS_CAPSELLA/SCRIPTS/")

# import data
df<-read.delim("../data/cluster_assignments.txt")
df<-as.data.frame(unique(df[,c("sample", "K", "assignment")]))

df1<-read.delim("../data/rad_data_samples_cw.txt")
names(df1)[1]<-"sample"

## merge datasets
df<-merge(df, df1, by="sample")

# extract kgcs
query<-na.omit(unique(df[,c("sample", "Lat", "Long")]))
names(query)<-c("sample", "Latitude", "Longitude")
clim<-as.data.frame(cbind(query, LookupCZ(query, res="coarse", rc=T)))
names(clim)[4]<-"zone"

df1<-merge(df, clim, by="sample")

##########
# K = 2
k2<-df1[df1$K==2,]

tbl<-as.data.frame(table(k2$assignment, k2$zone))
tbl$Var1<-factor(tbl$Var1, levels=c("Cluster1", "Cluster2", "Admixed"))
tbl$Var2<-as.character(tbl$Var2)
tbl<-tbl[!tbl$Var2=="Climate Zone info missing",]
tbl<-as.data.frame(cbind(tbl, str_split_fixed(tbl$Var2, "", 3)))
names(tbl)[4:6]<-c("1st", "2nd", "3rd")

# check by group
pal<-c("#D55E00", "#009E73", "#56B4E9", "#CC79A7", "dimgrey")
names(pal)<-c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Admixed")
s<-aggregate(Freq ~ Var1, data=tbl, FUN="sum")
names(s)<-c("Var1", "total")
tbl<-merge(tbl, s, by="Var1")
tbl$prop<-tbl$Freq/tbl$total
tbl<-tbl[tbl$prop > 0,]

p1<-ggplot(tbl, aes(x=Var2, y=prop, fill=Var1)) + 
  geom_bar(stat="identity", position="dodge") +
  theme_cowplot() +
  theme(legend.position="top") +
  scale_fill_manual(values=pal)
ggsave("~/Desktop/k2_summ.pdf", p1, height=4, width=8)

## major climate zone
tbl2<-aggregate(Freq ~ `1st` + Var1, data=tbl, FUN="sum")

p1<-ggplot(tbl2, aes(x=Var1, y=Freq, fill=`1st`)) + 
  geom_bar(stat="identity", position="fill") +
  theme_cowplot() +
  scale_fill_nejm()
ggsave("005A_k2_clim_group.pdf", p1, height=4, width=5)

## subplots per zone
tbl2<-aggregate(Freq ~ `2nd` + `1st` + Var1, data=tbl, FUN="sum")

p1<-ggplot(tbl2, aes(x=Var1, y=Freq, fill=`2nd`)) + 
  geom_bar(stat="identity", position="fill") +
  theme_cowplot() +
  scale_fill_nejm() +
  facet_grid(.~`1st`)
ggsave("005B_k2_clim_subgroup.pdf", p1, height=4, width=12)

##########
# K = 3
k2<-df1[df1$K==3,]

tbl<-as.data.frame(table(k2$assignment, k2$zone))
tbl$Var1<-factor(tbl$Var1, levels=c("Cluster1", "Cluster2", "Cluster3", "Admixed"))
tbl$Var2<-as.character(tbl$Var2)
tbl<-tbl[!tbl$Var2=="Climate Zone info missing",]
tbl<-as.data.frame(cbind(tbl, str_split_fixed(tbl$Var2, "", 3)))
names(tbl)[4:6]<-c("1st", "2nd", "3rd")

# check by group
s<-aggregate(Freq ~ Var1, data=tbl, FUN="sum")
names(s)<-c("Var1", "total")
tbl<-merge(tbl, s, by="Var1")
tbl$prop<-tbl$Freq/tbl$total
tbl<-tbl[tbl$prop > 0,]

p1<-ggplot(tbl, aes(x=Var2, y=prop, fill=Var1)) + 
  geom_bar(stat="identity", position="dodge") +
  theme_cowplot() +
  theme(legend.position="top") +
  scale_fill_manual(values=pal)
ggsave("~/Desktop/k3_summ.pdf", p1, height=4, width=8)

## major climate zone
tbl2<-aggregate(Freq ~ `1st` + Var1, data=tbl, FUN="sum")

p1<-ggplot(tbl2, aes(x=Var1, y=Freq, fill=`1st`)) + 
  geom_bar(stat="identity", position="fill") +
  theme_cowplot() +
  scale_fill_nejm()
ggsave("005C_k3_clim_group.pdf", p1, height=4, width=5)

## subplots per zone
tbl2<-aggregate(Freq ~ `2nd` + `1st` + Var1, data=tbl, FUN="sum")

p1<-ggplot(tbl2, aes(x=Var1, y=Freq, fill=`2nd`)) + 
  geom_bar(stat="identity", position="fill") +
  theme_cowplot() +
  scale_fill_nejm() +
  facet_grid(.~`1st`)
ggsave("005D_k3_clim_subgroup.pdf", p1, height=4, width=14)

# K = 4
k2<-df1[df1$K==4,]

tbl<-as.data.frame(table(k2$assignment, k2$zone))
tbl$Var1<-factor(tbl$Var1, levels=c("Cluster1", "Cluster2", "Cluster3", 
                                    "Cluster4", "Admixed"))
tbl$Var2<-as.character(tbl$Var2)
tbl<-tbl[!tbl$Var2=="Climate Zone info missing",]
tbl<-as.data.frame(cbind(tbl, str_split_fixed(tbl$Var2, "", 3)))
names(tbl)[4:6]<-c("1st", "2nd", "3rd")

# check by group
s<-aggregate(Freq ~ Var1, data=tbl, FUN="sum")
names(s)<-c("Var1", "total")
tbl<-merge(tbl, s, by="Var1")
tbl$prop<-tbl$Freq/tbl$total
tbl<-tbl[tbl$prop > 0,]

p1<-ggplot(tbl, aes(x=Var2, y=prop, fill=Var1)) + 
  geom_bar(stat="identity", position="dodge") +
  theme_cowplot() +
  theme(legend.position="top") +
  scale_fill_manual(values=pal)
ggsave("~/Desktop/k4_summ.pdf", p1, height=4, width=8)


## major climate zone
tbl2<-aggregate(Freq ~ `1st` + Var1, data=tbl, FUN="sum")

p1<-ggplot(tbl2, aes(x=Var1, y=Freq, fill=`1st`)) + 
  geom_bar(stat="identity", position="fill") +
  theme_cowplot() +
  scale_fill_nejm()
ggsave("005E_k4_clim_group.pdf", p1, height=4, width=6)

## subplots per zone
tbl2<-aggregate(Freq ~ `2nd` + `1st` + Var1, data=tbl, FUN="sum")

p1<-ggplot(tbl2, aes(x=Var1, y=Freq, fill=`2nd`)) + 
  geom_bar(stat="identity", position="fill") +
  theme_cowplot() +
  scale_fill_nejm() +
  facet_grid(.~`1st`)
ggsave("005F_k4_clim_subgroup.pdf", p1, height=3, width=16)

## quick climate map
library(pacman)
p_load(readxl, ggplot2, ggmap, maps, mapdata,viridis, kgc, 
       cowplot, magrittr, maptools, raster, rgeos, rasterVis, viridis, 
       scatterpie, tidyr, ggsci)

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

df1<-as.data.frame(cbind(df1, str_split_fixed(df1$zone, "", 3)))
p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df1, 
                  aes(x=Longitude.x, y=Latitude.x, color=`1`)) +
  coord_fixed(xlim=c(-20, 50), ylim=c(-40, 40), ratio=1) +
  theme_void() +
  scale_fill_npg() +
  theme(legend.position="top")
ggsave("~/Desktop/clim1_africa.pdf", p1, height=5.5, width=5)


p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df1, 
             aes(x=Longitude.x, y=Latitude.x, color=`2`)) +
  coord_fixed(xlim=c(-20, 50), ylim=c(-40, 40), ratio=1) +
  theme_void() +
  scale_fill_npg() +
  theme(legend.position="top")
ggsave("~/Desktop/clim2_africa.pdf", p1, height=5.5, width=5)

#df1<-df1[df1$Latitude.x >= -20 & df1$Latitude.x <=20,]
#df1<-df1[df1$Longitude.x >= 0 & df1$Longitude.x <=20,]
ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df1, 
             aes(x=Longitude.x, y=Latitude.x, color=Elevation)) +
  coord_fixed(xlim=c(-20, 20), ylim=c(0, 20), ratio=1) +
  theme_void() +
  scale_fill_viridis() +
  theme(legend.position="top")
