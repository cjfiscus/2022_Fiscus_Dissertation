#!/usr/bin/env Rscript
# plot het
# cjfiscus
# 2022-03-21

library(pacman)
p_load(ggplot2, cowplot, dplyr, reshape2)

setwd("~/Desktop/FIGURE_DRAFTS_COWPEA/")

# Heterozygosity plots
df<-read.delim("data/het.txt")

## plot He
p1<-ggplot(df, aes(x=He)) + geom_density() +
  theme_cowplot()
ggsave("015A_He.pdf", p1, height=4, width=4)

p1<-ggplot(df, aes(x=Ho)) + geom_density() +
  theme_cowplot()
ggsave("015B_Ho.pdf", p1, height=4, width=4)

m<-melt(df, id.vars="snp_id")

p1<-ggplot(m, aes(x=value, fill=variable)) + geom_histogram(color="black", alpha=0.5) + 
  theme_cowplot() +
  theme(legend.position="top")
ggsave("015C_He_Ho.pdf", p1, height=4, width=4)
