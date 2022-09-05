#!/usr/bin/env Rscript
# coverage filter
# cjfiscus

library(pacman)
p_load(data.table, ggplot2, matrixStats)

args = commandArgs(trailingOnly=TRUE)

# fun to import and format data 
importData<-function(path){
  counts<-fread(path)
  counts<-as.data.frame(counts)
  row.names(counts)<-counts$mer
  counts$mer<-NULL
  dim(counts)
  return(counts)
}

# import data
counts<-importData(args[1])

# tabulate sums and format table 
df<-as.data.frame(cbind(colnames(counts), colSums(counts)))
rm(counts) # purge
names(df)<-c("Library", "CountSum")
df$CountSum<-as.numeric(as.character(df$CountSum))

# calculate coverage
df$Coverage<-df$CountSum/150000000

# plot
g<-ggplot(df) + geom_density(aes(x=Coverage), fill="grey") + theme_classic() + 
	theme(text=element_text(size=16)) + xlab("Estimated Coverage")
ggsave("../results/filters/cov.png", g, height=5, width=5, units="in")

# threshold at 1X
df1<-df[df$Coverage < 1,]

write.table(df1$Library, args[3], sep="\t", quote=F, row.names=F, col.names=F)

# write out full table
write.table(df, args[2], sep="\t", quote=F, row.names=F)
