#!/usr/bin/env Rscript
library(pacman)
p_load(data.table)

# removes libraries to be filtered from K-mer count table

args = commandArgs(trailingOnly=TRUE)

# read in data 
counts<-fread(args[1])

# ids to rm
drop<-read.delim("../results/filters/filtered_similar_lines.txt")
drop$x<-as.character(drop$x)
drop<-drop$x

# filter
counts[,c(drop):=NULL]
dim(counts)

# write out
fwrite(counts, file=args[2], sep="\t", quote=F)
