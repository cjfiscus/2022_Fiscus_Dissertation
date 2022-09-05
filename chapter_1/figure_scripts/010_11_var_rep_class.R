#!/usr/bin/env Rscript
# Variability in repeat classes vs buscos
# cjfiscus
# 2021-12-09

library(pacman)
p_load(data.table, ggplot2, cowplot, matrixStats)

## read in data
df<-fread("~/Desktop/data/A_thal_kmer_abund_norm.txt")
df<-as.data.frame(df)
df1<-fread("~/Desktop/data/A_thal_kmer_abund_busco.txt")
df1<-as.data.frame(df1)
df1[,2:ncol(df1)]<-2^df1[,2:ncol(df1)]
rpt<-read.delim("~/Desktop/data/A_thaliana_Repeat_Database.txt")
rpt<-rpt[,c("ID", "CLASS", "SUBCLASS", "FAMILY")]

## subset data to repeats
df<-df[df$Feature %in% rpt$ID,]

## calc var ests per seq
### repeats
v<-as.data.frame(cbind(df$Feature, rowMedians(as.matrix(df[,2:ncol(df)])), 
                       rowMaxs(as.matrix(df[,2:ncol(df)])),
                       rowMins(as.matrix(df[,2:ncol(df)]))))
names(v)<-c("ID", "median", "max", "min")
v$range<-as.numeric(as.character(v$max))-as.numeric(as.character(v$min))
v$rom<-as.numeric(as.character(v$range))/as.numeric(as.character(v$median))

### buscos
v1<-as.data.frame(cbind(df1$Feature, rowMedians(as.matrix(df1[,2:ncol(df1)])), 
                       rowMaxs(as.matrix(df1[,2:ncol(df1)])),
                       rowMins(as.matrix(df1[,2:ncol(df1)]))))
names(v1)<-c("ID", "median", "max", "min")
v1$range<-as.numeric(as.character(v1$max))-as.numeric(as.character(v1$min))
v1$rom<-as.numeric(as.character(v1$range))/as.numeric(as.character(v1$median))

## merge data frames
v<-merge(v, rpt, by="ID")
v1$CLASS<-"BUSCO"
v1$SUBCLASS<-"BUSCO"
v1$FAMILY<-"BUSCO"
v<-as.data.frame(rbind(v, v1))

## produce plot
v<-v[!v$CLASS=="Other",]
v$CLASS<-factor(v$CLASS, levels=c("BUSCO", "Retrotransposon", "DNA transposon",
                                    "Satellite", "Simple Repeat"))
p1<-ggplot(v, aes(x=CLASS, y=log2(rom))) + geom_boxplot() +
  theme_cowplot() + ylab("log(sR)")
ggsave("~/Desktop/10_sr_by_class.pdf", p1, height=3, width=8)
##########
## breakout by class
sub<-v[!v$CLASS=="BUSCO",]

### retrotrans
ord<-as.data.frame(aggregate(rom ~ FAMILY + CLASS, data=sub, FUN=median))
ord<-ord[order(ord$rom),]
ord1<-ord[ord$CLASS=="Retrotransposon",]
sub1<-sub[sub$CLASS=="Retrotransposon",]
sub1$FAMILY<-factor(sub1$FAMILY, levels=ord1$FAMILY)
p1<-ggplot(sub1, aes(x=FAMILY, y=log2(rom))) + 
  geom_boxplot() +
  theme_cowplot() + ylab("log(sR)")
ggsave("~/Desktop/11A_retrotransposons.pdf", p1, height=3, width=8)

### DNA trans
ord1<-ord[ord$CLASS=="DNA transposon",]
sub1<-sub[sub$CLASS=="DNA transposon",]
sub1$FAMILY<-factor(sub1$FAMILY, levels=ord1$FAMILY)
p1<-ggplot(sub1, aes(x=FAMILY, y=log2(rom))) + 
  geom_boxplot() +
  theme_cowplot() + ylab("log(sR)")
ggsave("~/Desktop/11B_dnatransposons.pdf", p1, height=3, width=8)

### Satellites
sub1<-sub[sub$CLASS=="Satellite",]
ord<-as.data.frame(aggregate(rom ~ SUBCLASS, data=sub1, FUN=median))
ord<-ord[order(ord$rom),]
sub1$SUBCLASS<-factor(sub1$SUBCLASS, levels=ord$SUBCLASS)

p1<-ggplot(sub1, aes(x=SUBCLASS, y=log2(rom))) + geom_boxplot() +
  theme_cowplot() + ylab("log(sR)")
ggsave("~/Desktop/11C_satellites.pdf", p1, height=3, width=8)

### Simple Repeats
sub1<-sub[sub$CLASS=="Simple Repeat",]
ord<-as.data.frame(aggregate(rom ~ SUBCLASS, data=sub1, FUN=median))
ord<-ord[order(ord$rom),]
sub1$SUBCLASS<-factor(sub1$SUBCLASS, levels=ord$SUBCLASS)

p1<-ggplot(sub1, aes(x=SUBCLASS, y=log2(rom))) + geom_boxplot() +
  theme_cowplot() + ylab("log(sR)")
ggsave("~/Desktop/11D_simplerepeats.pdf", p1, height=3, width=8)
