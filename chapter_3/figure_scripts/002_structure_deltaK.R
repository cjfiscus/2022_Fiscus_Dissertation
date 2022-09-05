#!/usr/bin/env Rscript
# structure delta K plot
# cjfiscus
# 2022-03-04

library(pacman)
p_load(ggplot2, cowplot, dplyr)

setwd("~/Desktop/FIGURE_DRAFTS_COWPEA/SCRIPTS/")

# import data
df<-read.table("../DATA/logli.txt")
names(df)<-c("K", "logli")

## calc mean and sd of logli
d<-as.data.frame(aggregate(logli ~ K, data=df, FUN=mean))
names(d)<-c("K", "mlogli")
e<-as.data.frame(aggregate(logli ~ K, data=df, FUN=sd))
names(e)<-c("K", "slogli")
d<-merge(d, e, by="K")
d$delta<-d$mlogli-lag(d$mlogli)
d$delta2<-abs(lead(d$delta)-d$delta)
d$deltak<-d$delta2/d$slogli

### plot mean logli
p1<-ggplot(d, aes(x=K, y=mlogli)) + geom_line() + geom_point() +
  theme_cowplot() +
  geom_errorbar(aes(ymin=mlogli-slogli, ymax=mlogli+slogli), width=.2) +
  ylab("log likelihood") +
  scale_x_continuous(breaks=seq(1,10,1))
ggsave("002A_logli.pdf", p1, height=4, width=4)

### plot delta K
p1<-ggplot(d, aes(x=K, y=deltak)) + geom_line() + geom_point() +
  theme_cowplot() +
  scale_x_continuous(breaks=seq(1,10,1)) +
  ylab(expression(Delta*"K"))
ggsave("002B_admixture_deltaK.pdf", p1, height=4, width=4)
