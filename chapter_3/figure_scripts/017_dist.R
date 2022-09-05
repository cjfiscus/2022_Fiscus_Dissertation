#!/usr/bin/env Rscript
# plot dist
# cjfiscus
# 2022-03-21

library(pacman)
p_load(ggplot2, cowplot, dplyr, reshape2)

setwd("~/Desktop/FIGURE_DRAFTS_COWPEA/")

df<-read.table("data/plink.dist")
nm<-read.table("data/plink.dist.id")
names(df)<-nm$V1
row.names(df)<-nm$V2
xy <- t(combn(colnames(df), 2))
pw<-data.frame(xy, d=df[xy])

p1<-ggplot(pw, aes(x=d)) + geom_histogram(binwidth=1000) +
  theme_cowplot() +
  xlab("genetic distance")
ggsave("017_dist.pdf", p1, height=4, width=4)  
