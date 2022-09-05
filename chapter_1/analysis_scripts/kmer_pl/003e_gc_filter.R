#!/usr/bin/env Rscript
# extreme gc filter
# cjfiscus

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2)
args = commandArgs(trailingOnly=TRUE)

# read in data
df<-read.delim(args[1])

# define outlying values
low<-quantile(df$GC_content)[2] - 2*IQR(df$GC_content)
high<-quantile(df$GC_content)[4] + 2*IQR(df$GC_content)

# produce plot
g<-ggplot(df) + geom_density(aes(x=GC_content), fill="grey") + 
	geom_vline(aes(xintercept=low), alpha=0.5, color="red") + 
	geom_vline(aes(xintercept=high), alpha=0.5, color="red")
ggsave("../results/filters/gc_dist.png", g, height=5, width=5, units="in")

## apply filters 
df1<-df[df$GC_content < low,]
df2<-df[df$GC_content > high,]

df3<-rbind(df1, df2)

write.table(df3$Library, args[2], sep="\t", quote=F, row.names=F, col.names=F)
