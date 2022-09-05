#!/usr/bin/env Rscript
# folded af spectrum
# cjfiscus
# 2022-03-21

library(pacman)
p_load(ggplot2, cowplot, dplyr)

setwd("~/Desktop/FIGURE_DRAFTS_COWPEA/")

# allele freq
df<-read.table("data/plink.frq", header=T)

breaks<-seq(0,0.5, 0.1)
df$bin<-cut(df$MAF, breaks=breaks, include.lowest=F)
df<-na.omit(df)
p1<-ggplot(df, aes(x=bin)) + geom_bar() +
  theme_cowplot() + xlab("minor allele frequency") +
  ylab("Number of sites")
ggsave("014_maf.pdf", p1, height=4, width=4)
