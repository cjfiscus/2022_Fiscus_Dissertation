#!/usr/bin/env Rscript
# gc corr filter
# cjfiscus

args = commandArgs(trailingOnly=TRUE)

# new data
df<-read.table(args[1])

# threshold
h<-quantile(df$V2)[4] + 1.5*IQR(df$V2)
l<-quantile(df$V2)[2] - 1.5*IQR(df$V2)

# new bad list 
df2<-df[df$V2 > h,]
df3<-df[df$V2 < l,]

out<-rbind(df2, df3)

write.table(out$V1, args[2], sep="\t", quote=F, col.names=F, row.names=F)
