#!/usr/bin/env Rscript
# pipeline filter figs
# cjfiscus
# 2021-06-28

setwd("~/Desktop/FINAL_FIGS/")
library(pacman)
p_load(ggplot2, cowplot, ggdendro, patchwork)

# % Mapping filter
df<-read.table("./data/mapstats_genome.txt")
names(df)<-c("lib", "per")

p1<-ggplot(df, aes(x=per)) + geom_density() +
  theme_cowplot() + xlab("% reads mapped") + 
  geom_vline(aes(xintercept=90), color="red")
ggsave("filter1.pdf", p1, height=3, width=3)

# coverage filter
df<-read.delim("./data/coverage.txt")
df<-df[1:1257,]
p2<-ggplot(df, aes(x=Coverage)) + geom_density() +
  theme_cowplot() +
  xlab("estimated coverage") +
  geom_vline(aes(xintercept=1), color="red")
ggsave("filter2.pdf", p2, height=3, width=3)

# Extreme GC content
df<-read.delim("./data/gc_content_lib.txt")
h<-quantile(df$GC_content)[4] + 1.5*IQR(df$GC_content)
l<-quantile(df$GC_content)[2] - 1.5*IQR(df$GC_content)
p3<-ggplot(df, aes(x=GC_content)) + geom_density() +
  theme_cowplot() +
  xlab("GC content") +
  geom_vline(aes(xintercept=l), color="red") +
  geom_vline(aes(xintercept=h), color="red")
ggsave("filter3.pdf", p3, height=3, width=3)

# near identical lines
df<-read.table("./data/dist.dist.gz")
nm<-read.table("./data/dist.dist.id")
names(df)<-nm$V1
row.names(df)<-nm$V1

xy <- t(combn(colnames(df), 2))
pw<-data.frame(xy, dist=df[xy])

p4<-ggplot(pw, aes(x=dist)) + geom_density() +
  theme_cowplot() +
  xlab("genetic distance") +
  geom_vline(aes(xintercept=100000), color="red") +
  scale_x_continuous(labels=function(x)x/1000000) +
  #xlab("genetic distance (X 106)") +
  labs(x = expression ("genetic distance (X"~10^6~")"))
  
ggsave("filter4.pdf", p4, height=3, width=3)

# GC/count correlation
df<-read.table("./data/gc_cor_post.txt")
names(df)<-c("id", "corr")

h1<-quantile(df$corr)[4]+1.5*IQR(df$corr)
l1<-quantile(df$corr)[2]-1.5*IQR(df$corr)
scaleFUN <- function(x) sprintf("%.2f", x)
p5<-ggplot(df, aes(x=corr)) + geom_density() +
  xlab("Spearman's correlation") + 
  theme_cowplot() +
  geom_vline(aes(xintercept=l1), color="red") +
  geom_vline(aes(xintercept=h1), color="red") +
  scale_x_continuous(labels=scaleFUN)
ggsave("filter5.pdf", p5, height=3, width=3)

## arrange plots 
grid<-p1 + p2 + p3 + p4 + p5 + plot_layout(ncol=3) +
  plot_annotation(tag_levels="A")
ggsave("filters_all.pdf", grid, height=6, width=9)
