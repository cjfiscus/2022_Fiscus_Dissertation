#!/usr/bin/env Rscript
# pca of pop structure proportion variance plot
# cjfiscus
# 2022-03-07

library(pacman)
p_load(ggplot2, cowplot, readxl, kgc, ggsci, stringr)

setwd("~/Desktop/FIGURE_DRAFTS_COWPEA/SCRIPTS/")

# import data
## pc scores
df<-read.delim("../data/scores_pc.txt")
df$sample<-row.names(df)

## cowpea passport info
df1<-read_excel("../data/IITA-Cowpea collection.xls")
df1<-df1[,c("Accession name", "Latitude", "Longitude", "Biological status of accession")]
names(df1)[1]<-"sample"
query<-na.omit(df1[,c("sample", "Latitude", "Longitude")])
clim<-as.data.frame(cbind(query, LookupCZ(query, res="coarse", rc=T)))
names(clim)[4]<-"zone"

## cluster assignments
df2<-read.delim("../data/cluster_assignments.txt")
df2<-as.data.frame(unique(df2[,c("sample", "K", "assignment")]))

## merge data frames for plotting
df<-merge(df, df1, by="sample", all.x=T)
df<-merge(df, clim, by="sample", all.x=T)
df1<-merge(df, df2, by="sample", all.x=T)

## proportion var
prop_var<-read.delim("../data/prop_var_pc.txt")

# plots
## all points black A
p1<-ggplot(df, aes(x=Axis1, y=Axis2)) + geom_point() +
  theme_cowplot() +
  xlab(paste0("PC 1 (", round(prop_var[1,3], digits=3)*100, "%)")) +
  ylab(paste0("PC 2 (", round(prop_var[2,3], digits=3)*100, "%)"))

ggsave("011A_pca_pc1_pc2.pdf", p1, height=4, width=4)

## colored by germplasm type B
p1<-ggplot(df, aes(x=Axis1, y=Axis2, color=`Biological status of accession`)) + 
  geom_point() +
  theme_cowplot() +
  xlab(paste0("PC 1 (", round(prop_var[1,3], digits=3)*100, "%)")) +
  ylab(paste0("PC 2 (", round(prop_var[2,3], digits=3)*100, "%)")) +
  theme(legend.position="right")
ggsave("011B_pca_pc1_pc2_biostatus.pdf", p1, height=4, width=6)

## colored by climate zone major D
zones<-as.data.frame(cbind(df$sample, str_split_fixed(df$zone, "", 3)))
names(zones)<-c("sample", "zone_major", "zone_minor", "zone_minor2")
df<-merge(df, zones, by="sample", all.x=T)
df$zone_major<-factor(df$zone_major, levels=c("A", "B", "C"))

p1<-ggplot(df, aes(x=Axis1, y=Axis2, color=zone_major)) + 
  geom_point() +
  theme_cowplot() +
  xlab(paste0("PC 1 (", round(prop_var[1,3], digits=3)*100, "%)")) +
  ylab(paste0("PC 2 (", round(prop_var[2,3], digits=3)*100, "%)")) +
  theme(legend.position="right")  +
  scale_color_nejm()
ggsave("011C_pca_pc1_pc2_clim_group.pdf", p1, height=4, width=6)

## colored by climate zone minor E
df$zone_minor<-factor(df$zone_minor, levels=c("f", "m", "s", "S", "w", "W"))

p1<-ggplot(df, aes(x=Axis1, y=Axis2, color=zone_minor)) + 
  geom_point() +
  theme_cowplot() +
  xlab(paste0("PC 1 (", round(prop_var[1,3], digits=3)*100, "%)")) +
  ylab(paste0("PC 2 (", round(prop_var[2,3], digits=3)*100, "%)")) +
  theme(legend.position="right")  +
  scale_color_nejm()

ggsave("011D_pca_pc1_pc2_clim_subgroup.pdf", p1, height=4, width=6)

## colored by K = 2 D
df1$assignment<-factor(df1$assignment, levels=c("Cluster1", "Cluster2", "Cluster3", "Cluster4",
                                    "Admixed"))
k2<-df1[df1$K==2,]

p1<-ggplot() +
  geom_point(data=k2[k2$assignment=="Admixed",], aes(x=Axis1, y=Axis2), color="gray69") +
  geom_point(data=k2[!k2$assignment=="Admixed",], aes(x=Axis1, y=Axis2, color=assignment)) +
  geom_point() +
  theme_cowplot() +
  xlab(paste0("PC 1 (", round(prop_var[1,3], digits=3)*100, "%)")) +
  ylab(paste0("PC 2 (", round(prop_var[2,3], digits=3)*100, "%)")) +
  theme(legend.position="right")  +
  scale_color_npg()

ggsave("011E_pca_pc1_pc2_struc_k2.pdf", p1, height=4, width=6)

## colored by K = 3 E
k2<-df1[df1$K==3,]

p1<-ggplot() +
  geom_point(data=k2[k2$assignment=="Admixed",], aes(x=Axis1, y=Axis2), color="gray69") +
  geom_point(data=k2[!k2$assignment=="Admixed",], aes(x=Axis1, y=Axis2, color=assignment)) +
  geom_point() +
  theme_cowplot() +
  xlab(paste0("PC 1 (", round(prop_var[1,3], digits=3)*100, "%)")) +
  ylab(paste0("PC 2 (", round(prop_var[2,3], digits=3)*100, "%)")) +
  theme(legend.position="right")  +
  scale_color_npg()

ggsave("011F_pca_pc1_pc2_struc_k3.pdf", p1, height=4, width=6)

## colored by K = 4 F

k2<-df1[df1$K==4,]

p1<-ggplot() +
  geom_point(data=k2[k2$assignment=="Admixed",], aes(x=Axis1, y=Axis2), color="gray69") +
  geom_point(data=k2[!k2$assignment=="Admixed",], aes(x=Axis1, y=Axis2, color=assignment)) +
  geom_point() +
  theme_cowplot() +
  xlab(paste0("PC 1 (", round(prop_var[1,3], digits=3)*100, "%)")) +
  ylab(paste0("PC 2 (", round(prop_var[2,3], digits=3)*100, "%)")) +
  theme(legend.position="right")  +
  scale_color_npg()

ggsave("011G_pca_pc1_pc2_struc_k4.pdf", p1, height=4, width=6)
