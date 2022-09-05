#!/usr/bin/env Rscript
# Figure 1
# cjfiscus
# 2022-03-21

library(pacman)
p_load(ggplot2, cowplot, ggsci, tidyr, stringr, 
       readxl, patchwork, reshape2, plyr,
       readxl, ggplot2, ggmap, maps, mapdata,viridis, kgc, 
       cowplot, magrittr, maptools, raster, rgeos, rasterVis, viridis)

setwd("~/Desktop/MANUSCRIPTS/FIGURE_DRAFTS_COWPEA/SCRIPTS/")

# plot map of collection sites colored by bio status
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

df$`Biological status of accession`<-revalue(df$`Biological status of accession`, 
                                             c("Breeding / Research material"="Breeding",
                                               "Traditional cultivar / Landrace"="Landrace"))

## pal and shapes
pal<-c("#4477AA", "#228833", "#AA3377", "#EE6677")
names(pal)<-c("Breeding", "Landrace", "Weedy", "Wild")
s<-c(9, 16, 15, 17)
names(s)<-c("Breeding", "Landrace", "Weedy", "Wild")

map<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3)+
  geom_point(data=df[df$`Biological status of accession`=="Landrace",], 
             aes(x=Longitude, y=Latitude, color=`Biological status of accession`,
                 shape=`Biological status of accession`), 
             size=1.5, alpha=0.4) +
  geom_point(data=df[!df$`Biological status of accession`=="Landrace",], 
             aes(x=Longitude, y=Latitude, color=`Biological status of accession`,
                 shape=`Biological status of accession`), 
             size=1.55, alpha=0.8) +
  theme_void() +
  coord_fixed(xlim=c(-25, 55), ylim=c(-35, 35), ratio=1) +
  scale_color_manual(values=pal) +
  scale_shape_manual(values=s) +
  theme(legend.position="none")
ggsave("~/Desktop/map.pdf", map, height=6, width=6)

leg<-get_legend(map + theme(legend.position="right"))
ggsave("~/Desktop/mapleg.pdf", leg, height=4, width=4)
##########

# B folded allele frequency spectrum
df<-read.table("../DATA/af.frq.strat", header=T)

breaks<-seq(0,0.5, 0.05)
df$bin<-cut(df$MAF, breaks=breaks, include.lowest=F)
df<-na.omit(df)

df$CLST<-revalue(df$CLST, c("Traditional"="Landrace"))
pal<-c("#4477AA", "#228833", "#AA3377", "#EE6677")
names(pal)<-c("Breeding", "Landrace", "Weedy", "Wild")

p1<-ggplot(df, aes(x=bin, fill=CLST)) + geom_bar(color="black", position="dodge") +
  theme_cowplot() + xlab("minor allele frequency") +
  ylab("Number of sites") +
  scale_x_discrete(breaks=unlist(split(levels(df$bin),
                                       rep(1:4, length=length(levels(df$bin))))[1])) +
  scale_fill_manual(values=pal) +
  theme(legend.position=c(0.70,0.90), legend.direction="vertical",
              legend.title=element_blank())

# C het dist
df<-read.delim("../DATA/het.txt")

m<-melt(df, id.vars="snp_id")

pal<-c("#0078AD", "#00C4A0")
names(pal)<-c("Ho", "He")

#ggplot(df, aes(x=He)) + geom_histogram()

p2<-ggplot(m, aes(x=value, fill=variable)) + 
  geom_histogram(color="black",alpha=0.5, position="identity", binwidth = 0.025) + 
  theme_cowplot() +
  theme(legend.position=c(0.5,1), legend.direction="horizontal",
        legend.title=element_blank()) +
  scale_fill_manual(values=pal) +
  scale_x_continuous(limits=c(-.025, 0.55), breaks=seq(0,0.5, 0.25))

# D Finbreeding
df<-read.table("../DATA/fstat.het", header=T)

p3<-ggplot(df, aes(x=F)) + geom_histogram(binwidth=0.05, 
                                          fill="#104E8B", color="black") +
  theme_cowplot() +
  scale_x_continuous(breaks=seq(-0.5, 1, 0.5))

# E distance
df<-read.table("../DATA/dist.dist")
nm<-read.table("../DATA/dist.dist.id")
names(df)<-nm$V1
row.names(df)<-nm$V2
xy <- t(combn(colnames(df), 2))
pw<-data.frame(xy, d=df[xy])

p4<-ggplot(pw, aes(x=d)) + geom_histogram(binwidth=1000, fill="#104E8B", color="black") +
  theme_cowplot() +
  xlab("genetic distance")
ggsave("~/Desktop/dist.pdf", p4, width=4.5, height=4)

##########
# arrange fig
fig<-p1 + p2 + p3 + p4 + plot_layout(ncol=2) + plot_annotation(tag_levels = 'A')
ggsave("~/Desktop/fig1_DRAFT.pdf", fig, width=8, height=8)

#ggsave("~/Desktop/fig1_DRAFT2.pdf", map, width=8, height=8)
