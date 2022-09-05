#!/usr/bin/env Rscript
# norm cov ests
# cjfiscus
# 2022-02-14

library(pacman)
p_load(dplyr, tidyr)

args=commandArgs(trailingOnly = T)

id<-args[1]

# import data
df<-read.table(paste0("/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/coverage/ests/",id, "_repeats.txt"))
names(df)<-c("CHR", "START", "END", "ID", "SCORE")
df$LEN<-df$END-df$START

## filter to > 50 bp
db<-read.delim("/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/data/A_thaliana_Repeat_Database.txt")
db<-db[,c("ID", "LENGTH")]
df1<-merge(df, db, by="ID")
df1$prop<-df1$LEN/df1$LENGTH
df1<-df1[df1$prop >=0.5,]
#df1<-df[df$LEN >=50,]

## sum per annotation
df<-df %>% group_by(ID) %>% summarize(total=sum(SCORE))
df1<-df1 %>% group_by(ID) %>% summarize(total=sum(SCORE))

# normalize by busco
bus<-read.table(paste0("/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/coverage/ests/",id, "_busco.txt"))
names(bus)<-c("CHR", "START", "END", "BUSCO_ID", "SCORE")
med<-median(bus$SCORE)
df$NORM<-df$total/med
df1$NORM<-df1$total/med
df$sample<-id
df1$sample<-id

# write out
write.table(df, paste0("/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/coverage/ests/",id, "_mnorm.txt"), sep="\t", quote=F, row.names=F, col.names=F)
write.table(df1, paste0("/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/results/coverage/ests/",id, "_mnorm_half.txt"), sep="\t", quote=F, row.names=F, col.names=F)
