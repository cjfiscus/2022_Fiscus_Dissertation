#!/usr/bin/env Rscript
# structure barplots
# cjfiscus
# 2022-03-04

library(pacman)
p_load(readxl, ggplot2, cowplot, ggsci, tidyr, stringr)

setwd("~/Desktop/FIGURE_DRAFTS_COWPEA/SCRIPTS/")

# import data
df<-read.delim("../data/cluster_assignments.txt")
df$sample<-toupper(df$sample)

df1<-read_excel("../data/IITA-Cowpea collection.xls")
df1<-df1[,c("Accession name", "Latitude", "Longitude", "Country of origin", "Biological status of accession")]
df1$`Accession name`<-toupper(df1$`Accession name`)
names(df1)[1]<-"sample"

## merge datasets
df<-merge(df, df1, by="sample")
df$variable<-NULL
df$value<-NULL
df<-as.data.frame(unique(df))

##########
# K = 2
k2<-df[df$K==2,]
tbl<-as.data.frame(table(k2$assignment, k2$`Biological status of accession`))
tbl$Var1<-factor(tbl$Var1, levels=c("Cluster1", "Cluster2", "Admixed"))

p1<-ggplot(tbl, aes(x=Var1, y=Freq, fill=Var2)) + geom_bar(stat="identity") +
  theme_cowplot() +
  theme(legend.position="right")
ggsave("004A_structure_biostatus_k2.pdf", p1, height=4, width=8)

##########
# K = 3
k3<-df[df$K==3,]
tbl<-as.data.frame(table(k3$assignment, k3$`Biological status of accession`))
tbl$Var1<-factor(tbl$Var1, levels=c("Cluster1", "Cluster2", "Cluster3", "Admixed"))

p1<-ggplot(tbl, aes(x=Var1, y=Freq, fill=Var2)) + geom_bar(stat="identity") +
  theme_cowplot() +
  theme(legend.position="right")
ggsave("004B_structure_biostatus_k3.pdf", p1, height=4, width=8)


# K = 4
k4<-df[df$K==4,]
tbl<-as.data.frame(table(k4$assignment, k4$`Biological status of accession`))
tbl$Var1<-factor(tbl$Var1, levels=c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Admixed"))

p1<-ggplot(tbl, aes(x=Var1, y=Freq, fill=Var2)) + geom_bar(stat="identity") +
  theme_cowplot() +
  theme(legend.position="right")
ggsave("004C_structure_biostatus_k4.pdf", p1, height=4, width=8)
