#!/usr/bin/env Rscript
# pca of pop structure proportion variacne plot
# cjfiscus
# 2022-03-07

library(pacman)
p_load(ggplot2, cowplot)

setwd("~/Desktop/FIGURE_DRAFTS_COWPEA/SCRIPTS/")

# import data
df<-read.delim("../data/prop_var_pc.txt")
df<-df[1:10,]

p1<-ggplot(df, aes(x=pc, y=proportion_variance*100)) + 
  geom_bar(stat="identity") +
  theme_cowplot() +
  xlab("Principal component") +
  ylab("Percent of variance explained") +
  scale_x_continuous(breaks=seq(1,10,1)) +
  scale_y_continuous(breaks=seq(0, 10,2.5), limits=c(0, 10))
ggsave("010_pca_per_var.pdf", p1, height=4, width=4)
