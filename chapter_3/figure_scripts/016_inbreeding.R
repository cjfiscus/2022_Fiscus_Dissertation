#!/usr/bin/env Rscript
# plot F
# cjfiscus
# 2022-03-21

library(pacman)
p_load(ggplot2, cowplot, dplyr, reshape2)

setwd("~/Desktop/FIGURE_DRAFTS_COWPEA/")

# inbreeding stat
df<-read.table("data/plink.het", header=T)

p1<-ggplot(df, aes(x=F)) + geom_histogram(binwidth=0.1) +
  theme_cowplot()
ggsave("016_F_inbreeding.pdf", p1, height=4, width=4)
