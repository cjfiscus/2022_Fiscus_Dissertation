#!/usr/bin/env Rscript
# plot rep vs non-rep stats
# cjfiscus
# 2021-12-17

setwd("~/Desktop")
library(pacman)
p_load(ggplot2, cowplot)

# import data
df<-read.table("~/Desktop/MANUSCRIPTS/athal_data/stats.txt")
names(df)<-c("stat", "class", "freq", "N", "K")

df$len<-ifelse(df$class=="repetitive",  21182514, 98650419)

## plot unique
df1<-df[df$stat=="unique",]

p1<-ggplot(df1, aes(x=K, y=freq/N, fill=class)) + geom_bar(stat="identity", position="dodge") +
  scale_x_continuous(breaks=seq(5,20,1)) +
  theme_cowplot() + ylab("unique / total K-mers") +
  theme(legend.position="top") + scale_fill_brewer(palette="Paired")
ggsave("~/Desktop/2_unique.pdf", p1, height=3, width=6)

## plot median abun
df2<-df[df$stat=="median_abun",]
df2$value<-(df2$freq*100000000)/df2$len

p2<-ggplot(df2, aes(x=K, y=log10(value), fill=class)) + 
  geom_bar(stat="identity", position="dodge") +
  scale_x_continuous(breaks=seq(5,20,1)) +
  theme_cowplot() + ylab("log (median count per 100 Mbp)") +
  theme(legend.position="top") + scale_fill_brewer(palette="Paired")
ggsave("2_abundance.pdf", p2, height=3, width=6)

library(patchwork)
out<-p1/p2 + theme(legend.position="none")
ggsave("~/Desktop/out.pdf", out, height=6, width=6)
