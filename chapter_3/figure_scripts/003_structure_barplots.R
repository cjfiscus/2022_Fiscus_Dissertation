#!/usr/bin/env Rscript
# structure barplots
# cjfiscus
# 2022-03-04

library(pacman)
p_load(ggplot2, cowplot, ggsci, tidyr, stringr)

setwd("~/Desktop/FIGURE_DRAFTS_COWPEA/SCRIPTS/")

# import data
df<-read.delim("../data/cluster_assignments.txt")

##########
# K = 2
k2<-df[df$K==2,]

## produce sort order
ord<-k2[,c("sample", "variable", "value")]
ord<-spread(ord, variable, value)
maxval<-apply(ord[,2:ncol(ord)], 1, max)
matchval<-vector(length=nrow(ord))
for(j in 1:nrow(ord)) matchval[j] <- match(maxval[j], ord[j, ])
ord$maxval<-maxval
ord$matchval<-matchval
ord<-ord[with(ord, order(matchval, -maxval)),]

k2$sample<-factor(k2$sample, levels=ord$sample)

p1<-ggplot(k2, aes(x=sample, y=value, fill=variable)) + 
  geom_bar(stat="identity", position="stack") +
  theme_cowplot() +
  theme(legend.position="top",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Proportion of ancestry") +
  xlab("accession") +
  scale_fill_npg()

ggsave("003A_structure_barplot_k2.pdf", p1, height=3, width=10)

##########
# K = 3
k3<-df[df$K==3,]

## produce sort order
ord<-k3[,c("sample", "variable", "value")]
ord<-spread(ord, variable, value)
maxval<-apply(ord[,2:ncol(ord)], 1, max)
matchval<-vector(length=nrow(ord))
for(j in 1:nrow(ord)) matchval[j] <- match(maxval[j], ord[j, ])
ord$maxval<-maxval
ord$matchval<-matchval
ord<-ord[with(ord, order(matchval, -maxval)),]

k3$sample<-factor(k3$sample, levels=ord$sample)

### plot
p1<-ggplot(k3, aes(x=sample, y=value, fill=variable)) + 
  geom_bar(stat="identity", position="stack") +
  theme_cowplot() +
  theme(legend.position="top",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Proportion of ancestry") +
  xlab("accession") +
  scale_fill_npg()
ggsave("003B_structure_barplot_k3.pdf", p1, height=3, width=10)

# K = 4
k4<-df[df$K==4,]

## produce sort order
ord<-k4[,c("sample", "variable", "value")]
ord<-spread(ord, variable, value)
maxval<-apply(ord[,2:ncol(ord)], 1, max)
matchval<-vector(length=nrow(ord))
for(j in 1:nrow(ord)) matchval[j] <- match(maxval[j], ord[j, ])
ord$maxval<-maxval
ord$matchval<-matchval
ord<-ord[with(ord, order(matchval, -maxval)),]

k4$sample<-factor(k4$sample, levels=ord$sample)

### plot
p1<-ggplot(k4, aes(x=sample, y=value, fill=variable)) + 
  geom_bar(stat="identity", position="stack") +
  theme_cowplot() +
  theme(legend.position="top",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Proportion of ancestry") +
  xlab("accession") +
  scale_fill_npg()
ggsave("003C_structure_barplot_k4.pdf", p1, height=3, width=10)
