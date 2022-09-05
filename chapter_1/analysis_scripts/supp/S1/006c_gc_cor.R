#!/usr/bin/env Rscript
library(pacman)
p_load(data.table)

df<-read.delim("../results/12mers_annot.txt")
out<-as.data.frame(cor(df[,1], df[,3:4], method="spearman"))
write.table(out, "../results/gc_cor.txt", sep="\t", quote=F, row.names=F)
