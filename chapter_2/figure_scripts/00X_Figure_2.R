#!/usr/bin/env Rscript
# Figure 2
# cjfiscus
# 2022-07-18

library(pacman)
p_load(readxl, ggplot2, ggmap, maps, mapdata,viridis, kgc, 
       cowplot, magrittr, maptools, raster, rgeos, rasterVis, viridis, 
       scatterpie, tidyr, ggsci, ggsci, stringr, ggpubr, ggbeeswarm)

setwd("~/Desktop/MANUSCRIPTS/FIGURE_DRAFTS_CAPSELLA/SCRIPTS/")

# panel A world map
## import data
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
df2$assignment<-factor(df2$assignment, levels=c("Cluster1", "Cluster2", "Admixed"))

pal<-c("#4DBBD5FF", "#E64B35FF", "dimgrey")
names(pal)<-c("Cluster1", "Cluster2", "Admixed")

p1<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
                  aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(ratio=1) +
  theme_void() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        legend.position="none") +
  coord_cartesian(clip = "off") +
  scale_color_manual(values=pal) +
  geom_hline(aes(yintercept=42.35), linetype="dashed", alpha=0.5) +
  geom_hline(aes(yintercept=34.65), linetype="dashed", alpha=0.5) +
  geom_rug(data=df2[df2$assignment=="Cluster1",], 
           aes(y=Lat, color=assignment), outside = F, sides = "r", alpha=0.5, 
           size=0.1) +
  geom_rug(data=df2[df2$assignment=="Cluster2",], 
           aes(y=Lat, color=assignment), outside = TRUE, sides = "r", alpha=0.5, 
           size=0.1)

p1a<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat, color=assignment)) +
  coord_fixed(ratio=1) +
  theme_void() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        legend.position="none") +
  coord_cartesian(clip = "off") +
  scale_color_manual(values=pal) +
  geom_hline(aes(yintercept=42.35), linetype="dashed", alpha=0.5) +
  geom_hline(aes(yintercept=34.65), linetype="dashed", alpha=0.5) +
  geom_rug(data=df2[df2$assignment=="Cluster1",], 
           aes(y=Lat, color=assignment), outside = F, sides = "r", alpha=0.5, 
           size=0.1) +
  geom_rug(data=df2[df2$assignment=="Admixed",], 
           aes(y=Lat, color=assignment), outside = TRUE, sides = "r", alpha=0.5, 
           size=0.1)

ggsave("~/Desktop/A2.pdf", p1a, height=6, width=8)

## plain plot for presentation
p<-ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), 
               fill='NA', color='black', size=0.3) +
  geom_point(data=df2, 
             aes(x=Long, y=Lat), color="chocolate2") +
  coord_fixed(ratio=1) +
  theme_void() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        legend.position="none") +
  coord_cartesian(clip = "off") +
  scale_color_manual(values=pal)
ggsave("~/Desktop/samplemap.pdf", p, height=6, width=8)

##########
## B PCA
df<-read.delim("../data/scores_pc.txt")
df$sample<-row.names(df)
df<-merge(df, df2, by="sample")
pervar<-read.delim("../data/prop_var_pc.txt")

p2<-ggplot(df, aes(x=Axis1, y=Axis2, color=assignment)) + geom_point() +
  theme_cowplot() + theme(legend.position="none") +
  scale_color_manual(values=pal) + 
  xlab(paste0("PC 1 (", round(pervar[1,3], digits=2)*100, "%)")) +
  ylab(paste0("PC 2 (", round(pervar[2,3], digits=2)*100, "%)"))
##########
## C clim dists
df<-read.delim("../data/kgc_zones.txt")
df<-merge(df, df2, by="sample")
df<-as.data.frame(cbind(df, str_split_fixed(df$zone, "", 3)))
names(df)[8]<-"group"
tbl<-as.data.frame(table(df$assignment, df$group))
tbl<-tbl[tbl$Freq > 0,]
names(tbl)<-c("assignment", "group", "Freq")
s<-read.delim("../data/kcc.txt")
tbl<-merge(tbl, s, by="group")

pal2<-c("#2CA02CFF", "#8C564BFF", "#FF7F0EFF", "#9467BDFF", "#1F77B4FF")
names(pal2)<-c("tropical", "dry", "temperate","continental", "polar")
tbl$assignment<-factor(tbl$assignment, levels=c("Cluster1", "Admixed", "Cluster2"))

p3<-ggplot(tbl, aes(x=assignment, y=Freq, fill=zone2)) + 
  geom_bar(stat="identity", position=position_dodge2(preserve="single")) +
  theme_cowplot() + theme(legend.position="top", axis.title.x=element_blank()) +
  scale_fill_manual(values=pal2)
##########
## D-X common garden phenos
df<-df1[df1$K==2,]
df$variable<-NULL
df$value<-NULL
df<-as.data.frame(unique(df))
df<-df[,c("sample", "assignment", 
          "FlowStart.day_after_sowing._Fam_mean",
          "number_of_basal_branches_fam_mean",
          "plant_height_in_cm_mean_fam",
          "germination_percentage",
          "seed_weight_in_mg_per_50_seeds",
          "Relative_Genome_size_to_petroselinum_crispum",
          "LeafType")]
names(df)<-c("sample", "assignment",
             "Days to flowering",
             "Number of basal branches",
             "Plant height (cm)",
             "Proportion germinated",
             "Weight / 50_seeds (mg)",
             "Genome size (pg)",
             "Leaf type")

### common garden phenotypes
df$assignment<-factor(df$assignment, levels=c("Cluster1", "Admixed", "Cluster2"))
df$`Weight / 50_seeds (mg)`<-as.numeric(df$`Weight / 50_seeds (mg)`)
df$`Genome size (pg)`<-ifelse(df$`Genome size (pg)` < 0.18, NA, df$`Genome size (pg)`)
j<-0

for (i in unique(names(df)[3:8])){
  print(i)
  sub<-df[,c("sample", "assignment", i)]
  sub<-na.omit(sub)
  cnt<-as.data.frame(table(sub$assignment))
  names(cnt)<-c("assignment", "n")
  sub<-merge(sub, cnt, by="assignment")
  comps<-list(c(as.character(cnt[1,2]), as.character(cnt[3,2])))
  sub$n<-factor(sub$n, levels=cnt$n)
  
  p<-ggplot(sub, aes(x=n, y=as.numeric(sub[,3]), color=assignment)) + 
    geom_quasirandom(varwidth = TRUE) +
    stat_summary(fun=median,geom="point", color="black", size=3, shape=15, 
                 fill="white") +
    stat_compare_means(comparisons=comps, label = "p.signif") +
    scale_color_manual(values=pal) + theme_cowplot() + ylab(as.character(i)) +
    theme(legend.position="none",
          axis.title.x= element_blank(),
          axis.ticks.x=element_blank()) +
    scale_y_continuous( breaks=pretty_breaks())
  ggsave(paste0("~/Desktop/",j,".pdf"), p, height=3, width=2)
  assign(paste0("b",j), p)
  j<-j+1
 }

#### leaftype
sub<-df[,c("sample", "assignment", "Leaf type")]
sub<-na.omit(sub)
sub$`Leaf type`<-ifelse(sub$`Leaf type`=="unsicher", "unknown", sub$`Leaf type`)
tbl<-as.data.frame(table(sub$assignment, sub$`Leaf type`))
names(tbl)<-c("assignment", "type", "Freq")

pal3<-c("#D62728FF", "#E377C2FF", "#17BECFFF", "#BCBD22FF", "#7F7F7FFF")
names(pal3)<-c("het", "rho", "sim", "ten", "unknown")

p4<-ggplot(tbl, aes(x=assignment, y=Freq, fill=type)) + 
  geom_bar(stat="identity", position=position_dodge2(preserve="single")) +
  theme_cowplot() + 
  theme(legend.position="top") +
  scale_fill_manual(values=pal3) +
  theme(axis.title.x = element_blank())

##########
# export plots
ggsave("~/Desktop/A.pdf", p1, height=6, width=8)
ggsave("~/Desktop/B.pdf", p2 + theme(legend.position="none"), height=3, width=3)
ggsave("~/Desktop/C.pdf", p3 + theme(legend.position="none"), height=3, width=3)
ggsave("~/Desktop/D.pdf", p4 + theme(legend.position="none"), height=3, width=3)

leg<-get_legend(p2+ theme(legend.position="right"))
ggsave("~/Desktop/leg1.pdf", leg, height=4, width=6)
leg<-get_legend(p3 + theme(legend.position="right"))
ggsave("~/Desktop/leg2.pdf", leg, height=4, width=6)
leg<-get_legend(p4 + theme(legend.position="right"))
ggsave("~/Desktop/leg3.pdf", leg, height=4, width=6)
