#!/usr/bin/env Rscript
# calc spearman corr between kmer GC and counts
# cjfiscus

library(pacman)
p_load(data.table, stringr)

args = commandArgs(trailingOnly=TRUE)
# args[1] is input file, args[2] is output file, args[3] is array id

# SLURM_TASK_ARRAY_ID -> col to read in
array<-as.integer(args[3])

# import kmer and column
df<-fread(args[1], select=c(1, array))

# calculate gc
df$gc<-str_count(df$mer, pattern="G|C")

# calc cor
out<-as.data.frame(cbind(names(df)[2], cor(df[,2], df[,3], method="spearman")))

# write out
write.table(out, args[2], sep="\t", col.names=F, row.names=F, quote=F)
