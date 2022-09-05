#!/usr/bin/env Rscript
# filter near identical lines
# cjfiscus
# 2021-02-21

library(pacman)
p_load(ggplot2, cowplot, ggdendro)

args = commandArgs(trailingOnly=TRUE)
# args[1] is output prefix

# import data
df<-read.table("dist.dist.gz")
nms<-read.table("dist.dist.id")
names(df)<-nms$V1
row.names(df)<-nms$V1

# restrict to lines
keep<-read.table(paste0(args[1], "keep_id.txt"))

df<-df[row.names(df)%in%keep$V1, names(df)%in%keep$V1]

# plot distribution
xy <- t(combn(colnames(df), 2))
pw<-data.frame(xy, dist=df[xy])

p1<-ggplot(pw, aes(x=dist)) + geom_density() +
  theme_cowplot() +
  geom_vline(aes(xintercept=100000))
ggsave(paste0(args[1], "dist.jpeg"), p1, height=4, width=4)

# cluster lines
hc<-hclust(as.dist(df))

p1<-ggdendrogram(hc) +
  geom_hline(aes(yintercept=100000), color="red", linetype="dotted")
ggsave(paste0(args[1], "dendro.jpeg"), p1, height=8, width=8)

# statictreecut
grp<-as.data.frame(cutree(hc, h=100000))
names(grp)<-"group"
grp$id<-row.names(grp)
write.table(grp, paste0(args[1], "groups.txt"), sep="\t", quote=F, row.names=F)

## check number in each group
tbl<-as.data.frame(table(grp$group))

## plot frequency
p1<-ggplot(tbl, aes(x=Freq)) + geom_bar() +
  theme_cowplot()
ggsave(paste0(args[1], "groups.jpeg"), p1, height=4, width=4)

idlst<-data.frame()
# sample one from each cluster
for (i in unique(grp$group)){
  sub<-grp[grp$group==i,]
  sub<-sub[sample(nrow(sub), 1),]
  idlst<-as.data.frame(rbind(idlst, sub))
}

# write out near identical lines
out<-grp[!grp$id %in% idlst$id,]

write.table(out$id, paste0(args[1], "filtered_similar_lines.txt"), quote=F, row.names=F)
