#!/usr/bin/env Rscript
# calculate R2 of columns 2 and 3
# cjfiscus
args<-commandArgs(trailingOnly=T)
options(stringsAsFactors = F)
library(pacman)
p_load(data.table)

df<-fread(args[1])
K<-args[2]

corr<-cbind(K, cor(log10(df[,2]+1), log10(df[,3]+1), use="pairwise.complete.obs")^2)
write.table(corr, args[3], quote=F, row.names=F, col.names=F, sep="\t")
