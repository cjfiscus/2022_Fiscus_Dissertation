#!/usr/bin/env Rscript
library(pacman)
p_load(data.table, tidyr)
args = commandArgs(trailingOnly=TRUE)

nm<-args[1]
df<-fread(nm)

## calculate spearman correlation and format table
values<-as.data.frame(sapply(df[,2:ncol(df)], cor, df$TAIR10, method="spearman"))
names(values)<-c("correlation")
values$sample<-row.names(values)

# parse sample name into multiple columns 
values<-values %>% separate(sample, c("coverage", "sample"), sep="_")

# remove TAIR10 corr with itself
values<-na.omit(values)

mer<-gsub("mers","", unlist(strsplit(nm, "[./]"))[6])

# add K to table 
values$K<-mer

# write out
out_name<-paste0("../results/tables/", mer,"mers_cor.txt")
write.table(values, out_name, sep="\t", row.names=F, quote=F, col.names=F)
